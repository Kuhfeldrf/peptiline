from django.urls import include, path
from django.contrib import admin

from . import views

urlpatterns = [
    path("health/", views.health_check, name="health_check"),
    path("admin/", admin.site.urls),
    # PeptiLine landing page is both the site root and the historical
    # /peptiline/ path -- see docs/SPLIT_PLAN.md section 6 (URL continuity):
    # the hosted MBPDB app redirects /peptiline/ here, so the path must keep
    # resolving even though it's also the standalone site's home page.
    path("", views.peptiline_landing, name="peptiline_landing"),
    path("peptiline/", views.peptiline_landing),
    path("supplementals/", views.peptiline_supplementals, name="peptiline_supplementals"),
    path("about_us/", views.about_us, name="about_us"),
    path("data_transformation/", include("data_transformation.urls")),
    path("data_analysis/", include("data_analysis.urls", namespace="data_analysis")),
    path("heatmap/", include("heatmap_viz.urls", namespace="heatmap_viz")),
]
