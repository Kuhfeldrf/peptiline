"""
Single-process startup bootstrap: migrate + clearsessions + optional
superuser creation. Mirrors the MBPDB monolith's peptide/management/commands
/bootstrap.py so start.sh only pays one interpreter start + django.setup()
on container boot.
"""
import os

from django.contrib.auth import get_user_model
from django.core.management import call_command
from django.core.management.base import BaseCommand


class Command(BaseCommand):
    help = "Run startup migrations, clear expired sessions, and ensure the configured superuser exists."

    def handle(self, *args, **options):
        self.stdout.write("Running database migrations...")
        try:
            call_command("migrate", "--run-syncdb", interactive=False)
        except Exception as exc:
            self.stderr.write(f"WARNING: Migrations failed: {exc}")

        call_command("clearsessions")

        username = os.environ.get("DJANGO_SUPERUSER_USERNAME")
        password = os.environ.get("DJANGO_SUPERUSER_PASSWORD")
        if username and password:
            User = get_user_model()
            if not User.objects.filter(username=username).exists():
                self.stdout.write("Creating superuser...")
                User.objects.create_superuser(
                    username, os.environ.get("DJANGO_SUPERUSER_EMAIL", ""), password
                )
                self.stdout.write("Superuser created successfully")
            else:
                self.stdout.write("Superuser already exists")
