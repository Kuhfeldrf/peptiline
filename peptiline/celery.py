from __future__ import absolute_import, unicode_literals
import os
from celery import Celery
from celery.signals import worker_init

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "peptiline.settings")

app = Celery("peptiline", broker="redis://127.0.0.1:6379/0")
app.config_from_object("django.conf:settings", namespace="CELERY")


@worker_init.connect
def worker_init(**kwargs):
    try:
        import pwd

        celery_user = pwd.getpwnam("celery_user")
        os.setuid(celery_user.pw_uid)
    except (ImportError, KeyError):
        pass


app.autodiscover_tasks()
