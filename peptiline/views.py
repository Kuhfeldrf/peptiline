from django.core.cache import cache
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render


def check_progress(request, task_id):
    """Poll a Celery task's progress, as written by run_blast_search_task /
    fetch_uniprot_task (data_transformation/tasks.py) into the cache keys
    below. Ported from MBPDB's peptide/views.py::check_progress, trimmed to
    the fields those two tasks actually populate.
    """
    progress = cache.get(f'progress_{task_id}', 0)
    size = cache.get(f'size_{task_id}')
    elapsed_time = cache.get(f'elapsed_time_{task_id}', 0.0)
    status = cache.get(f'status_{task_id}', 'in_progress')

    if status == 'in_progress':
        if progress > 0 and size:
            percent_progress = progress / size * 100
            estimated_time_remaining = (elapsed_time * 100 / percent_progress) - elapsed_time
        else:
            percent_progress = 0.0
            estimated_time_remaining = None

        return JsonResponse({
            'task_id': task_id,
            'percent_progress': percent_progress,
            'progress': progress,
            'size': size,
            'elapsed_time': elapsed_time,
            'estimated_time_remaining': estimated_time_remaining,
            'status': status,
        })

    response_data = {
        'task_id': task_id,
        'status': status,
        'percent_progress': 100 if status == 'complete' else 0.0,
        'progress': progress,
        'size': size,
        'elapsed_time': elapsed_time,
    }
    if status == 'failed':
        response_data['error'] = cache.get(f'error_{task_id}', '')
    return JsonResponse(response_data)


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
