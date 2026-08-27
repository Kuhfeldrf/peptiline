"""
Django forms for the data transformation wizard.
"""
from django import forms


class PeptidomicUploadForm(forms.Form):
    """Step 1: Upload peptidomic data file (single or merge mode)."""
    peptidomic_file = forms.FileField(
        label='Peptidomic Data File',
        help_text='CSV, TSV, or XLSX file containing peptidomic results',
        required=False,
        widget=forms.ClearableFileInput(attrs={
            'accept': '.csv,.tsv,.txt,.xlsx,.xls',
            'class': 'form-control-file',
        }),
    )
    functional_file = forms.FileField(
        label='Functional Annotation Data File (Optional)',
        help_text='CSV, TSV, or XLSX file with pre-computed functional annotation results',
        required=False,
        widget=forms.ClearableFileInput(attrs={
            'accept': '.csv,.tsv,.txt,.xlsx,.xls',
            'class': 'form-control-file',
        }),
    )
    similarity_threshold = forms.IntegerField(
        label='Functional Similarity Threshold',
        initial=80,
        required=False,
        min_value=0,
        max_value=100,
        widget=forms.Select(
            choices=[(i, f'{i}%') for i in range(0, 101, 10)],
            attrs={'class': 'form-control', 'style': 'width: 200px;'},
        ),
    )

    def clean(self):
        cleaned_data = super().clean()
        has_single = bool(cleaned_data.get('peptidomic_file'))
        has_merge = bool(self.files.getlist('merge_files'))
        if not has_single and not has_merge:
            raise forms.ValidationError(
                'Please upload a peptidomic data file or multiple files to merge.'
            )
        return cleaned_data


class FastaUploadForm(forms.Form):
    """Step 3: Upload a custom FASTA file to supplement the protein database."""
    fasta_file = forms.FileField(
        label='Custom FASTA Protein Sequences',
        help_text='FASTA file with protein sequences for mapping',
        widget=forms.ClearableFileInput(attrs={
            'accept': '.fasta,.fa,.faa,.txt',
            'class': 'form-control-file',
        }),
    )


class GroupUploadForm(forms.Form):
    """Step 2: Upload group definition JSON."""
    group_file = forms.FileField(
        label='Group Definition File',
        help_text='JSON file mapping group names to column lists',
        widget=forms.ClearableFileInput(attrs={
            'accept': '.json',
            'class': 'form-control-file',
        }),
    )
