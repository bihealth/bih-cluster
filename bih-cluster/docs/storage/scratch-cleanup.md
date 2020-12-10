# Automated Cleanup of Scratch Volumes

The `scratch` volumes are automatically cleaned up nightly with the following mechanism.

- The `scratch` volume of each user, group, and project is crawled to identifies files that are older than two weeks.
- These files are moved into a directory that is named by the current date `YYYY-MM-DD` into the same sub folders as it was below `scratch`.
- These scratch directories are removed after two more weeks.

The cleanup code is a (probably updated) version of the following code:

```python
#!/usr/bin/env python3
"""Tool to cleanup the scratch directories on BIH HPC."""

import argparse
import datetime
import logging
import os
from os.path import join, dirname, normpath, islink, exists, isdir
import pathlib
import re
import shutil
import sys

#: The maximum age of files (by mtime) before they are put into the thrash bin.
MAX_AGE_SCRATCH = datetime.timedelta(days=14)

#: The maximum age of trash cans to keep around.
MAX_AGE_TRASH = datetime.timedelta(days=14)

#: Name of trash can directory.
TRASH_DIR_NAME = "BIH_TRASH"


def setup_trash(args, now):
    """Setup the trash directory."""
    trash_path = join(args.scratch_dir, TRASH_DIR_NAME, now.strftime("%Y-%m-%d"))
    logging.debug("Today's trash dir is %s", trash_path)
    logging.debug("  - mkdir -p %s", trash_path)
    logging.debug("  - chown root:root %s", dirname(trash_path))
    logging.debug("  - chmod 0755 %s", dirname(trash_path))
    logging.debug("  - chown root:root %s", trash_path)
    logging.debug("  - chmod 0755 %s", trash_path)
    if not args.dry_run:
        # Ensure that today's trash dir either does not exist or is a directory.
        if exists(dirname(trash_path)) and not isdir(dirname(trash_path)):
            os.unlink(dirname(trash_path))
        elif exists(trash_path) and not isdir(trash_path):
            os.unlink(dirname(trash_path))
        os.makedirs(trash_path, mode=0o775, exist_ok=True)
        os.chown(dirname(trash_path), 0, 0)
        os.chmod(dirname(trash_path), 0o755)
        os.chown(trash_path, 0, 0)
        os.chmod(trash_path, 0o755)
    return trash_path

def move_files(args, now, trash_path):
    """Move files into the trash directory."""
    logging.debug("  Walking %s recursively", args.scratch_dir)
    for root, dirs, files in os.walk(args.scratch_dir):
        # Do not go into trash directory itself.
        if root == args.scratch_dir:
            dirs[:] = [d for d in dirs if d != TRASH_DIR_NAME]
        # Only create output directory once.
        dir_created = False
        # Walk files.
        for name in files:
            path = join(root, name)
            logging.debug("    considering %s", path)
            age = now - datetime.datetime.fromtimestamp(os.lstat(path).st_mtime)
            if age > MAX_AGE_SCRATCH:
                local = path[len(args.scratch_dir) + 1 :]
                dest = join(trash_path, local)
                logging.debug("    - chown root:root %s", path)
                logging.debug("    - chmod 0644 %s", path)
                if not args.dry_run:
                    os.lchown(path, 0, 0)
                    if not islink(path):
                        os.chmod(path, 0o644)
                if not dir_created:
                    dir_created = True
                    logging.debug("    - mkdir -p %s", dirname(dest))
                    if not args.dry_run:
                        os.makedirs(dirname(dest), mode=0o775, exist_ok=True)
                logging.debug("    - %s -> %s", path, dest)
                if not args.dry_run:
                    os.rename(path, dest)


def empty_trash(args, now):
    """Empty the trash cans."""

    def log_error(_func, path, exc_info):
        logging.warn("could not delete path %s: %s", path, exc_info)

    logging.debug("  Emptying the trash cans...")
    trash_base = join(args.scratch_dir, TRASH_DIR_NAME)
    if os.path.exists(trash_base):
        entries = os.listdir(trash_base)
    else:
        entries = []
    for entry in entries:
        if re.match(r"^\d\d\d\d-\d\d-\d\d", entry):
            path = join(args.scratch_dir, TRASH_DIR_NAME, entry)
            logging.info("    considering %s", path)
            folder_date = datetime.datetime.strptime(entry, "%Y-%m-%d")
            if now - folder_date > MAX_AGE_TRASH:
                logging.info("    - rm -rf %s", path)
                if not args.dry_run:
                    shutil.rmtree(path, onerror=log_error)
    logging.debug("  ... done emptying the trash cans.")


def run(args):
    """Actually execute the scratch directory cleanup."""
    now = datetime.datetime.now()  # use one start time
    logging.info("Starting to cleanup %s...", args.scratch_dir)
    trash_path = setup_trash(args, now)
    move_files(args, now, trash_path)
    empty_trash(args, now)
    logging.info("... done with cleanup.")


def main(argv=None):
    """Main entry point into the program."""
    parser = argparse.ArgumentParser()
    parser.add_argument("scratch_dir", metavar="SCRATCH_DIR", help="Path to the scratch directory")
    parser.add_argument("--verbose", "-v", dest="verbosity", default=1, action="count")
    parser.add_argument("--quiet", "-q", dest="quiet", default=False, action="store_true")
    parser.add_argument(
        "--dry-run",
        "-n",
        dest="dry_run",
        default=False,
        action="store_true",
        help="Do not actually perform the actions",
    )

    if not shutil.rmtree.avoids_symlink_attacks:
        raise RuntimeError("Cannot execute with rmtree on unsafe platforms!")

    args = parser.parse_args(argv)
    args.scratch_dir = normpath(args.scratch_dir)

    logging.basicConfig(format="%(asctime)-15s %(levelname)3.3s %(message)s")
    logger = logging.getLogger()
    if args.quiet:
        logger.setLevel(logging.WARNING)
    elif args.verbosity == 1:
        logger.setLevel(logging.INFO)
    elif args.verbosity > 1:
        logger.setLevel(logging.DEBUG)

    return run(args)


if __name__ == "__main__":
    sys.exit(main())
```