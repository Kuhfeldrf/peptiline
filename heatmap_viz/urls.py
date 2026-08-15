from django.urls import path
from . import views

app_name = 'heatmap_viz'

urlpatterns = [
    path('', views.index, name='index'),
    path('upload/', views.upload, name='upload'),
    path('transfer-from-dt/', views.transfer_from_dt, name='transfer_from_dt'),
    path('transfer-from-da/', views.transfer_from_da, name='transfer_from_da'),
    path('fetch-sequence/', views.fetch_sequence, name='fetch_sequence'),
    path('apply-fasta/', views.apply_fasta, name='apply_fasta'),
    path('get-specific-options/', views.get_specific_options, name='get_specific_options'),
    path('signal-peptide-info/', views.signal_peptide_info, name='signal_peptide_info'),
    path('plot/', views.plot, name='plot'),
    path('download-plot/', views.download_plot, name='download_plot'),
    path('download-static/', views.download_static_image, name='download_static_image'),
    path('save-compact-plot/', views.save_compact_plot, name='save_compact_plot'),
    path('download-compact-plot/', views.download_compact_plot, name='download_compact_plot'),
]
