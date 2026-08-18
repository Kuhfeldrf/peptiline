"""
UniProt Integration Module

This module provides functions and classes to integrate UniProt functionality 
into the peptide visualization notebooks.
"""

import re
import pandas as pd
from IPython.display import display, HTML
import ipywidgets as widgets
from .uniprot_client import UniProtClient

class UniProtIntegration:
    """
    Class to add UniProt integration to existing peptide visualization notebooks.
    
    This class can be used to add UniProt search functionality to either:
    1. Data_Visualization.ipynb
    2. Heatmap_Visualization_widget_volia.ipynb
    
    Usage:
        # Import the module
        from utils.uniprot_integration import UniProtIntegration
        
        # Initialize with a reference to your data transformation object
        uniprot_integration = UniProtIntegration(data_transformer)
        
        # Add the UI components
        uniprot_integration.add_uniprot_ui()
        
        # Whenever you want to check for missing proteins:
        uniprot_integration.check_missing_proteins()
    """
    
    def __init__(self, data_transformer):
        """
        Initialize the UniProt integration.
        
        Args:
            data_transformer: The data transformation object from the notebook
        """
        self.data_transformer = data_transformer
        self.uniprot_client = None
        self.uniprot_search_enabled = False
        self.missing_proteins = set()
        self.output_area = None
        
        # Initialize UI components
        self.uniprot_search = None
        self.uniprot_id_input = None
        self.uniprot_search_button = None
        self.batch_fetch_button = None
    
    def add_uniprot_ui(self, output_area=None):
        """
        Add UniProt UI components to the notebook.
        
        Args:
            output_area: Optional output area for displaying messages
        """
        if output_area is not None:
            self.output_area = output_area
        else:
            # Try to use the output area from the data transformer
            self.output_area = getattr(self.data_transformer, 'output_area', None)
            
            # Create a new output area if none exists
            if self.output_area is None:
                self.output_area = widgets.Output()
        
        # Create UI components
        self.uniprot_search = widgets.Checkbox(
            value=False,
            description='Use UniProt for missing proteins',
            layout=widgets.Layout(width='300px'),
            style={'description_width': 'initial'}
        )
        
        self.uniprot_id_input = widgets.Text(
            value='',
            placeholder='Enter UniProt ID (e.g., P01308)',
            description='UniProt ID:',
            disabled=True,
            layout=widgets.Layout(width='300px'),
            style={'description_width': 'initial'}
        )
        
        self.uniprot_search_button = widgets.Button(
            description='Search UniProt',
            disabled=True,
            button_style='primary',
            tooltip='Search UniProt for protein sequence',
            layout=widgets.Layout(width='150px')
        )
        
        self.batch_fetch_button = widgets.Button(
            description='Fetch All Missing',
            disabled=True,
            button_style='info',
            tooltip='Fetch all missing proteins from UniProt',
            layout=widgets.Layout(width='150px')
        )
        
        # Create a UI section for UniProt integration
        uniprot_box = widgets.VBox([
            widgets.HTML("<h3><u>Import from UniProt:</u></h3>"),
            self.uniprot_search,
            widgets.HBox([
                self.uniprot_id_input,
                self.uniprot_search_button,
                self.batch_fetch_button
            ])
        ], layout=widgets.Layout(margin='20px 0'))
        
        # Register event handlers
        self.uniprot_search.observe(self._on_uniprot_search_change, names='value')
        self.uniprot_search_button.on_click(self._on_uniprot_search_click)
        self.batch_fetch_button.on_click(self._on_batch_fetch_click)
        
        # Display the UI components
        display(uniprot_box)
        display(self.output_area)
    
    def _on_uniprot_search_change(self, change):
        """Handle UniProt search checkbox toggle"""
        self.uniprot_search_enabled = change['new']
        
        if change['new']:
            self.uniprot_id_input.disabled = False
            self.uniprot_search_button.disabled = False
            self.batch_fetch_button.disabled = False
            
            # Initialize the UniProt client
            if self.uniprot_client is None:
                try:
                    self.uniprot_client = UniProtClient()
                    with self.output_area:
                        display(HTML('<b style="color:green;">UniProt client initialized successfully.</b>'))
                        
                    # Check for missing proteins if we already have merged data
                    self.check_missing_proteins()
                except Exception as e:
                    with self.output_area:
                        display(HTML(f'<b style="color:red;">Error initializing UniProt client: {str(e)}</b>'))
                    self.uniprot_search.value = False
        else:
            self.uniprot_id_input.disabled = True
            self.uniprot_search_button.disabled = True
            self.batch_fetch_button.disabled = True
    
    def _on_uniprot_search_click(self, b):
        """Handle UniProt search button click"""
        if not self.uniprot_client:
            try:
                self.uniprot_client = UniProtClient()
            except Exception as e:
                with self.output_area:
                    display(HTML(f'<b style="color:red;">Error: UniProt client could not be initialized: {str(e)}</b>'))
                return
                
        protein_id = self.uniprot_id_input.value.strip()
        if not protein_id:
            with self.output_area:
                display(HTML('<b style="color:red;">Please enter a UniProt ID.</b>'))
            return
            
        with self.output_area:
            self.output_area.clear_output()
            display(HTML(f'<b>Searching UniProt for ID: {protein_id}...</b>'))
            
            try:
                name, species, sequence = self.uniprot_client.fetch_protein_info_with_sequence(protein_id)
                
                if sequence:
                    # Get the appropriate dictionary based on the notebook
                    if hasattr(self.data_transformer, 'proteins_dic'):
                        protein_dict = self.data_transformer.proteins_dic
                    elif hasattr(self.data_transformer, 'protein_dict'):
                        protein_dict = self.data_transformer.protein_dict
                    else:
                        with self.output_area:
                            display(HTML('<b style="color:red;">Error: Cannot find protein dictionary in data transformer.</b>'))
                        return
                    
                    # Store the protein data
                    protein_dict[protein_id] = {
                        "name": name if name else protein_id,
                        "sequence": sequence,
                        "species": species if species else "Unknown"
                    }
                    
                    display(HTML(
                        f'<b style="color:green;">Successfully imported protein from UniProt:</b><br>' +
                        f'<b>ID:</b> {protein_id}<br>' +
                        f'<b>Name:</b> {name if name else "Unknown"}<br>' +
                        f'<b>Species:</b> {species if species else "Unknown"}<br>' +
                        f'<b>Sequence length:</b> {len(sequence)} amino acids'
                    ))
                    
                    # Remove from missing proteins if it was there
                    if protein_id in self.missing_proteins:
                        self.missing_proteins.remove(protein_id)
                    
                    # Clear input field
                    self.uniprot_id_input.value = ''
                else:
                    display(HTML(f'<b style="color:red;">No sequence found for UniProt ID: {protein_id}</b>'))
            except Exception as e:
                display(HTML(f'<b style="color:red;">Error fetching data from UniProt: {str(e)}</b>'))
    
    def _on_batch_fetch_click(self, b):
        """Handle batch fetch button click"""
        self.fetch_batch_from_uniprot()
    
    def fetch_batch_from_uniprot(self):
        """Fetch a batch of missing proteins from UniProt"""
        if not self.missing_proteins or not self.uniprot_client:
            with self.output_area:
                display(HTML('<b style="color:orange;">No missing proteins to fetch.</b>'))
            return
            
        with self.output_area:
            self.output_area.clear_output()
            display(HTML(f'<b>Fetching batch of {len(self.missing_proteins)} proteins from UniProt...</b>'))
            
            try:
                # Get the appropriate dictionary based on the notebook
                if hasattr(self.data_transformer, 'proteins_dic'):
                    protein_dict = self.data_transformer.proteins_dic
                elif hasattr(self.data_transformer, 'protein_dict'):
                    protein_dict = self.data_transformer.protein_dict
                else:
                    display(HTML('<b style="color:red;">Error: Cannot find protein dictionary in data transformer.</b>'))
                    return
                
                # Use the batch fetch method
                protein_ids = list(self.missing_proteins)
                results = self.uniprot_client.fetch_proteins_batch(protein_ids)
                
                success_count = 0
                for protein_id in protein_ids:
                    try:
                        if protein_id in results:
                            name, species = results[protein_id]
                            
                            # We need to fetch the sequence separately
                            _, _, sequence = self.uniprot_client.fetch_protein_info_with_sequence(protein_id)
                            
                            if sequence:
                                protein_dict[protein_id] = {
                                    "name": name if name else protein_id,
                                    "sequence": sequence,
                                    "species": species if species else "Unknown"
                                }
                                self.missing_proteins.remove(protein_id)
                                success_count += 1
                                display(HTML(f'<span style="color:green;">Successfully fetched {protein_id}: {name}</span>'))
                            else:
                                display(HTML(f'<span style="color:orange;">No sequence found for {protein_id}</span>'))
                        else:
                            display(HTML(f'<span style="color:orange;">No data found for {protein_id}</span>'))
                    except Exception as e:
                        display(HTML(f'<span style="color:orange;">Error processing {protein_id}: {str(e)}</span>'))
                
                if success_count > 0:
                    display(HTML(f'<b style="color:green;">Successfully fetched {success_count} out of {len(protein_ids)} missing proteins.</b>'))
                else:
                    display(HTML(f'<b style="color:red;">Failed to fetch any proteins. Please try individual search.</b>'))
            except Exception as e:
                display(HTML(f'<b style="color:red;">Error in batch fetch: {str(e)}</b>'))
    
    def check_missing_proteins(self):
        """Check for proteins in merged data that are missing from the protein dictionary"""
        # Skip if UniProt search is not enabled
        if not self.uniprot_search_enabled:
            return
        
        # Get the merged dataframe
        if hasattr(self.data_transformer, 'merged_df'):
            merged_df = self.data_transformer.merged_df
        else:
            with self.output_area:
                display(HTML('<b style="color:red;">Error: Cannot find merged dataframe in data transformer.</b>'))
            return
        
        if merged_df is None or merged_df.empty:
            return
        
        if 'Master Protein Accessions' not in merged_df.columns:
            return
        
        # Get the appropriate dictionary based on the notebook
        if hasattr(self.data_transformer, 'proteins_dic'):
            protein_dict = self.data_transformer.proteins_dic
        elif hasattr(self.data_transformer, 'protein_dict'):
            protein_dict = self.data_transformer.protein_dict
        else:
            with self.output_area:
                display(HTML('<b style="color:red;">Error: Cannot find protein dictionary in data transformer.</b>'))
            return
        
        # Get unique protein IDs from merged data
        unique_protein_ids = set(merged_df['Master Protein Accessions'].dropna().unique())
        
        # Find missing proteins
        self.missing_proteins = unique_protein_ids - set(protein_dict.keys())
        
        if not self.missing_proteins:
            with self.output_area:
                display(HTML('<b style="color:green;">All proteins in the data have sequence information.</b>'))
            return
        
        with self.output_area:
            self.output_area.clear_output()
            missing_list = ", ".join(list(self.missing_proteins)[:10])
            if len(self.missing_proteins) > 10:
                missing_list += f" and {len(self.missing_proteins) - 10} more"
                
            display(HTML(
                f'<b style="color:orange;">Found {len(self.missing_proteins)} proteins missing sequence information.</b><br>' +
                f'Missing proteins: {missing_list}<br>' +
                f'<b>Use the buttons above to fetch these proteins from UniProt.</b>'
            ))


def add_uniprot_integration_to_dataframe(df, uniprot_client=None):
    """
    Add protein information (name, species) to a dataframe by querying UniProt.
    
    Args:
        df: pandas DataFrame with a 'Master Protein Accessions' column
        uniprot_client: Optional UniProtClient instance
        
    Returns:
        DataFrame with added 'protein_name' and 'protein_species' columns
    """
    if uniprot_client is None:
        uniprot_client = UniProtClient()
    
    # Check if required columns exist
    if 'Master Protein Accessions' not in df.columns:
        raise ValueError("DataFrame must have a 'Master Protein Accessions' column")
    
    # Add protein_name and protein_species columns if they don't exist
    if 'protein_name' not in df.columns:
        df['protein_name'] = ''
    if 'protein_species' not in df.columns:
        df['protein_species'] = ''
    
    # Get unique protein IDs
    unique_proteins = df['Master Protein Accessions'].dropna().unique()
    
    # Fetch protein information in batch
    if len(unique_proteins) > 0:
        try:
            results = uniprot_client.fetch_proteins_batch(list(unique_proteins))
            
            # Update dataframe with results
            for protein_id, (name, species) in results.items():
                mask = df['Master Protein Accessions'] == protein_id
                df.loc[mask, 'protein_name'] = name if name else protein_id
                df.loc[mask, 'protein_species'] = species if species else "Unknown"
        except Exception as e:
            print(f"Error fetching protein information: {str(e)}")
    
    return df 