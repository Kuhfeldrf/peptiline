from django.http import HttpResponse


class HealthCheckMiddleware:
    """
    Answer /health/ before CommonMiddleware's ALLOWED_HOSTS check runs.

    Azure Container Apps' startup/readiness/liveness probes hit the
    container directly over its internal cluster IP rather than the public
    FQDN, so the request's Host header is something like '100.100.2.79' --
    never a value that belongs in ALLOWED_HOSTS. Must be the first entry in
    MIDDLEWARE so it runs before SecurityMiddleware/CommonMiddleware touch
    request.get_host().
    """

    def __init__(self, get_response):
        self.get_response = get_response

    def __call__(self, request):
        if request.path == "/health/":
            return HttpResponse("ok", content_type="text/plain")
        return self.get_response(request)
