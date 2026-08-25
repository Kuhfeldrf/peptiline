from whitenoise.storage import CompressedManifestStaticFilesStorage


class LenientManifestStaticFilesStorage(CompressedManifestStaticFilesStorage):
    """
    Same as WhiteNoise's CompressedManifestStaticFilesStorage, except a
    {% static %} reference to a file that doesn't exist on disk (not just
    missing from the manifest -- manifest_strict=False alone doesn't cover
    this, since hashed_name() itself raises when it can't find the source
    file to hash) degrades to a broken link/image instead of raising
    ValueError and 500ing the whole page. This repo's vendoring pass (see
    docs/SPLIT_PLAN.md) didn't bring over ~50 marketing/screenshot/demo/
    supplemental assets referenced by base.html/peptiline_landing.html/
    peptiline_supplementals.html; sourcing them is separate follow-up work,
    not a reason every page should 500 in production until it's done.
    """

    manifest_strict = False

    def stored_name(self, name):
        try:
            return super().stored_name(name)
        except ValueError:
            return self.clean_name(name)
