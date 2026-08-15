from django.urls import path
from . import views

app_name = 'data_analysis'

urlpatterns = [
    path('', views.index, name='index'),
    path('upload/', views.upload, name='upload'),
    path('transfer-from-dt/', views.transfer_from_dt, name='transfer_from_dt'),
    path('plot/', views.plot, name='plot'),
    path('save-plot/', views.save_plot, name='save_plot'),
    path('download-plot/', views.download_plot, name='download_plot'),
    path('download-static/', views.download_static_image, name='download_static_image'),
]
