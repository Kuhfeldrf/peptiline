"""
UniProt Client Module

This module provides functions to query protein information from UniProt.
It includes functions to fetch single protein data or batch protein data.
"""

import requests
import time
from xml.etree import ElementTree


def fetch_uniprot_info(protein_id, fetch_sequence=False):
    """
    Fetch protein information from UniProt, prioritizing common names.
    
    Args:
        protein_id (str): UniProt accession ID.
        fetch_sequence (bool): Whether to fetch protein sequence. Defaults to False.
        
    Returns:
        tuple: (protein_common_name, species_common_name, sequence) or (None, None, None) if not found.
              If fetch_sequence is False, sequence will be None.
    """
    try:
        sequence = None
        # Try REST API first
        rest_url = f'https://rest.uniprot.org/uniprotkb/{protein_id}.json'
        response = requests.get(rest_url)
        
        if response.status_code == 200:
            data = response.json()
            
            # Get protein common name
            try:
                # Look for protein names
                names = data['proteinDescription']
                protein_name = None
                
                # Try to find a common/short name
                if 'shortNames' in names.get('recommendedName', {}):
                    protein_name = names['recommendedName']['shortNames'][0]['value']
                elif 'shortNames' in names.get('alternativeNames', [{}])[0]:
                    protein_name = names['alternativeNames'][0]['shortNames'][0]['value']
                else:
                    # Fallback to full name
                    protein_name = names.get('recommendedName', {}).get('fullName', {}).get('value')
                
                # Get species common name - improved logic
                organism_data = data.get('organism', {})
                species = None
                
                # Try all possible name types in order of preference
                name_types = ['common', 'scientific', 'synonym']
                for name_type in name_types:
                    for name in organism_data.get('names', []):
                        if name.get('type') == name_type:
                            species = name.get('value')
                            if species:  # If we found a valid name, break
                                break
                    if species:  # If we found a valid name, break
                        break
                
                # If still no species name, try lineage
                if not species and 'lineage' in organism_data:
                    species = organism_data['lineage'][-1]  # Get the most specific taxon
                
                # Get sequence if requested
                if fetch_sequence and 'sequence' in data:
                    sequence = data['sequence'].get('value')
                
                return protein_name, species, sequence
                
            except KeyError:
                pass  # Fall through to XML approach
        
        # Fall back to XML API
        xml_url = f'https://www.uniprot.org/uniprot/{protein_id}.xml'
        response = requests.get(xml_url)
        
        if response.status_code != 200:
            return None, None, None
            
        root = ElementTree.fromstring(response.content)
        ns = {'up': 'http://uniprot.org/uniprot'}
        
        # Get protein common name with fallbacks
        protein_name = None
        # Try short name first
        name_element = (
            root.find('.//up:recommendedName/up:shortName', ns) or
            root.find('.//up:alternativeName/up:shortName', ns) or
            root.find('.//up:recommendedName/up:fullName', ns) or
            root.find('.//up:submittedName/up:fullName', ns)
        )
        protein_name = name_element.text if name_element is not None else None
        
        # Get species name with improved logic
        species = None
        # Try different name types in order of preference
        name_types = ['common', 'scientific', 'synonym']
        for name_type in name_types:
            organism = root.find(f'.//up:organism/up:name[@type="{name_type}"]', ns)
            if organism is not None:
                species = organism.text
                if species:  # If we found a valid name, break
                    break
        
        # If still no species name, try lineage
        if not species:
            lineage = root.find('.//up:organism/up:lineage', ns)
            if lineage is not None:
                species = lineage.text.split(',')[-1].strip()  # Get the most specific taxon
        
        # Get sequence if requested
        if fetch_sequence:
            seq_element = root.find('.//up:sequence', ns)
            if seq_element is not None:
                sequence = seq_element.text
        
        if protein_name:
            # Clean up protein name - remove any "precursor" or similar suffixes
            protein_name = protein_name.split(' precursor')[0].split(' (')[0]
        else:
            return None, None, None
            
        return protein_name, species, sequence
        
    except Exception as e:
        print(f"Error fetching UniProt data for {protein_id}: {str(e)}")
        return None, None, None

def fetch_uniprot_info_batch(protein_ids, max_retries=3, timeout=30):
    """
    Fetch protein information for multiple proteins at once using UniProt's batch API.
    
    Args:
        protein_ids (list): List of UniProt accession IDs.
        max_retries (int, optional): Maximum number of retry attempts. Defaults to 3.
        timeout (int, optional): Request timeout in seconds. Defaults to 30.
        
    Returns:
        dict: Dictionary mapping protein IDs to (name, species) tuples.
    """
    if not protein_ids:
        return {}
    
    results = {}
    
    try:
        # Use the batch REST API endpoint
        batch_url = 'https://rest.uniprot.org/uniprotkb/search'
        
        # Create a query with OR conditions for each protein ID
        query = ' OR '.join([f'accession:{pid}' for pid in protein_ids])

        # Parameters for the request
        params = {
            'query': query,
            'format': 'json',
            'fields': 'accession,protein_name,organism_name',
            'size': len(protein_ids)  # Request all results in one response
        }
        
        # Make the request with timeout and retry logic
        retry_count = 0
        while retry_count < max_retries:
            try:
                response = requests.get(batch_url, params=params, timeout=timeout)
                
                # Handle rate limiting
                if response.status_code == 429:  # Too Many Requests
                    retry_after = int(response.headers.get('Retry-After', 5))
                    print(f"Rate limited. Waiting {retry_after} seconds...")
                    time.sleep(retry_after + 1)  # Add 1 second buffer
                    retry_count += 1
                    continue
                
                # Break if successful
                if response.status_code == 200:
                    break
                
                # Handle other errors
                response.raise_for_status()
                
            except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
                retry_count += 1
                wait_time = 2 ** retry_count  # Exponential backoff
                print(f"Request failed: {str(e)}. Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
        
        # Process the response
        if response.status_code == 200:
            data = response.json()
            
            for entry in data.get('results', []):
                accession = entry.get('primaryAccession')
                
                # Extract protein name
                protein_name = None
                protein_data = entry.get('proteinDescription', {})
                #{'recommendedName': {'fullName': {'value': 'Alpha-S2-casein'}}, 'contains': [{'recommendedName': {'fullName': {'value': 'Casocidin-1'}}, 'alternativeNames': [{'fullName': {'value': 'Casocidin-I'}}]}], 'flag': 'Precursor'}

                # Try to find a common/short name first
                if 'recommendedName' in protein_data:
                    if 'shortNames' in protein_data['recommendedName']:
                        protein_name = protein_data['recommendedName']['shortNames'][0]['value']
                    else:
                        protein_name = protein_data['recommendedName'].get('fullName', {}).get('value')
                
                # Try alternative names if no recommended name found
                if not protein_name and 'alternativeNames' in protein_data and protein_data['alternativeNames']:
                    if 'shortNames' in protein_data['alternativeNames'][0]:
                        protein_name = protein_data['alternativeNames'][0]['shortNames'][0]['value']
                    else:
                        protein_name = protein_data['alternativeNames'][0].get('fullName', {}).get('value')
                # Extract species name
                species = None
                organism_data = entry.get('organism')
                #{'scientificName': 'Bos taurus', 'commonName': 'Bovine', 'taxonId': 9913, 'lineage': ['Eukaryota', 'Metazoa', 'Chordata', 'Craniata', 'Vertebrata', 'Euteleostomi', 'Mammalia', 'Eutheria', 'Laurasiatheria', 'Artiodactyla', 'Ruminantia', 'Pecora', 'Bovidae', 'Bovinae', 'Bos']}
                if 'commonName' in str(organism_data):
                    protein_species = organism_data.get('commonName') 
                else:
                    protein_species = organism_data.get('scientificName') 

                # Clean up protein name if found
                if protein_name:
                    protein_name = protein_name.split(' precursor')[0].split(' (')[0]
                
                # Store the results
                if accession and protein_name and protein_species:
                    results[accession] = (protein_name, protein_species)# or "Unknown")
        
        return results
    
    except Exception as e:
        print(f"Error in batch fetch: {str(e)}")
        return {}

class UniProtClient:
    """
    A client class for interacting with the UniProt API.
    
    This class provides methods for fetching protein information from UniProt,
    with caching to avoid redundant requests.
    """
    
    def __init__(self):
        """Initialize the UniProt client with an empty cache."""
        self.cache = {}
        self.sequence_cache = {}
    
    
    def fetch_protein_info_with_sequence(self, protein_id):
        """
        Fetch information for a single protein including sequence, using cache if available.
        
        Args:
            protein_id (str): UniProt accession ID.
            
        Returns:
            tuple: (protein_name, species, sequence)
        """
        # Check sequence cache first
        if protein_id in self.sequence_cache:
            return self.sequence_cache[protein_id]
        
        # Check if we have basic info but not sequence
        if protein_id in self.cache:
            name, species = self.cache[protein_id]
            # Need to fetch sequence separately
            _, _, sequence = fetch_uniprot_info(protein_id, fetch_sequence=True)
            result = (name, species, sequence)
            self.sequence_cache[protein_id] = result
            return result
        
        # Fetch everything from UniProt
        name, species, sequence = fetch_uniprot_info(protein_id, fetch_sequence=True)
        
        # Cache the results
        self.cache[protein_id] = (name, species)
        self.sequence_cache[protein_id] = (name, species, sequence)
        
        return name, species, sequence
    
    def fetch_proteins_batch(self, protein_ids, max_retries=3, timeout=30):
        """
        Fetch information for multiple proteins in batch, using cache when possible.
        
        Args:
            protein_ids (list): List of UniProt accession IDs.
            max_retries (int, optional): Maximum number of retry attempts. Defaults to 3.
            timeout (int, optional): Request timeout in seconds. Defaults to 30.
            
        Returns:
            dict: Dictionary mapping protein IDs to (name, species) tuples.
        """
        # Filter out proteins that are already in the cache
        uncached_proteins = [pid for pid in protein_ids if pid not in self.cache]
        
        # Prepare results dictionary with cached values
        results = {pid: self.cache[pid] for pid in protein_ids if pid in self.cache}
        
        # If there are proteins not in the cache, fetch them
        if uncached_proteins:
            batch_results = fetch_uniprot_info_batch(uncached_proteins, max_retries, timeout)
            
            # Update cache and results
            for pid, info in batch_results.items():
                self.cache[pid] = info
                results[pid] = info
        
        return results 