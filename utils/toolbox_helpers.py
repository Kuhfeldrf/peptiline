"""
Standalone upload/cleanup helpers vendored from the parent MBPDB application's
peptide/toolbox.py. Only the two functions PeptiLine's own modules import are
included here; toolbox.py's other functions depend on MBPDB's proprietary
models and are not part of this package.
"""
import os
import shutil


def handle_uploaded_file(request_file, path):
    with open(path, 'wb') as destination:
        for chunk in request_file.chunks():
            destination.write(chunk)


def clear_temp_directory(directory_path):
    dirs = [f for f in os.scandir(directory_path) if f.is_dir()]
    dirs.sort(key=lambda x: os.path.getmtime(x.path), reverse=True)
    for dir_entry in dirs[25:]:
        try:
            shutil.rmtree(dir_entry.path)
        except Exception:
            pass
