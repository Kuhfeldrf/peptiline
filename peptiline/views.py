from django.http import HttpResponse
from django.shortcuts import render


def health_check(request):
    # Deliberately cheap (no DB hit) -- Container Apps probes hit this
    # through nginx -> gunicorn so a ready replica actually means Django
    # is serving, not just that nginx's socket is open.
    return HttpResponse("ok", content_type="text/plain")


def peptiline_landing(request):
    return render(request, "peptide/peptiline_landing.html")


def peptiline_supplementals(request):
    return render(request, "peptide/peptiline_supplementals.html")


def about_us(request):
    return render(request, "peptide/about_us.html")
