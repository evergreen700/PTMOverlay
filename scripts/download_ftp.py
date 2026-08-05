"""Download a file, a glob match, or a whole directory tree from an FTP/FTPS server.

Usage: download_ftp.py [--unzip GLOB] credentials.yaml remote_path local_path
"""

import argparse
import fnmatch
import ftplib
import os
import posixpath
import shutil
import ssl
import time
import yaml
import zipfile

TIMEOUT = 120
RETRIES = 4
GLOB_CHARS = "*?["


def connect(cred):
    """Log in, retrying: the server rate-limits new connections."""
    for attempt in range(1, RETRIES + 1):
        try:
            ftp = ftplib.FTP_TLS(timeout=TIMEOUT)
            ftp.connect(cred["address"], cred.get("port", 21))
            ftp.login(
                user=cred.get("username", "anonymous"),
                passwd=cred.get("password", ""),
            )
            # MassIVE (and most FTPS servers) refuse unencrypted data connections
            ftp.prot_p()
            # so that SIZE is answered instead of "550 not allowed in ASCII mode"
            ftp.voidcmd("TYPE I")
            return ftp
        except ftplib.all_errors as e:  # includes OSError/ssl.SSLError
            print(f"  connection attempt {attempt}/{RETRIES} failed: {e}", flush=True)
            if attempt == RETRIES:
                raise
            time.sleep(10 * attempt)


class Remote:
    def __init__(self, cred):
        self.cred = cred
        self.ftp = connect(cred)

    def reconnect(self):
        try:
            self.ftp.close()
        except Exception:
            pass
        self.ftp = connect(self.cred)

    def size(self, path):
        """File size in bytes, or None if path is a directory/missing."""
        # listings leave the connection in ASCII mode, where SIZE is refused
        self.ftp.voidcmd("TYPE I")
        try:
            return self.ftp.size(path)
        except ftplib.error_perm:
            return None

    def listdir(self, path):
        """(name, is_dir) pairs for a remote directory."""
        try:
            return [
                (name, facts.get("type") == "dir")
                for name, facts in self.ftp.mlsd(path)
                if name not in (".", "..")
            ]
        except ftplib.error_perm:
            pass  # server without MLSD support
        entries = []
        for entry in self.ftp.nlst(path):
            name = posixpath.basename(entry)
            if name in (".", ".."):
                continue
            entries.append((name, self.size(posixpath.join(path, name)) is None))
        return entries

    def resolve(self, path):
        """Expand a glob in the last path component to a single remote path."""
        base = posixpath.basename(path)
        if not any(c in base for c in GLOB_CHARS):
            return path
        parent = posixpath.dirname(path)
        matches = [n for n, _ in self.listdir(parent) if fnmatch.fnmatch(n, base)]
        if len(matches) != 1:
            raise SystemExit(
                f"ERROR: {path} matched {len(matches)} entries ({', '.join(sorted(matches))}); "
                "expected exactly 1"
            )
        return posixpath.join(parent, matches[0])

    def download(self, remote, local):
        for attempt in range(1, RETRIES + 1):
            try:
                expected = self.size(remote)
                with open(local, "wb") as outFile:
                    self.ftp.retrbinary(f"RETR {remote}", outFile.write)
                got = os.path.getsize(local)
                if expected is not None and got != expected:
                    raise OSError(f"got {got} of {expected} bytes")
                return
            except (ftplib.error_temp, ftplib.error_proto, ssl.SSLError, OSError, EOFError) as e:
                print(f"  attempt {attempt}/{RETRIES} failed: {e}", flush=True)
                if attempt == RETRIES:
                    raise
                time.sleep(5 * attempt)
                self.reconnect()

    def download_tree(self, remote, local):
        """Download a directory recursively; returns the local paths written."""
        os.makedirs(local, exist_ok=True)
        written = []
        for name, is_dir in sorted(self.listdir(remote)):
            remote_entry = posixpath.join(remote, name)
            local_entry = os.path.join(local, name)
            if is_dir:
                written += self.download_tree(remote_entry, local_entry)
            else:
                print(f"Downloading {remote_entry} -> {local_entry}", flush=True)
                self.download(remote_entry, local_entry)
                written.append(local_entry)
        return written


def extract_zips(paths, pattern):
    """Flatten matching members out of each downloaded .zip, then drop the archive."""
    for archive in paths:
        if not archive.endswith(".zip"):
            continue
        with zipfile.ZipFile(archive) as z:
            members = [
                m for m in z.namelist()
                if fnmatch.fnmatch(posixpath.basename(m), pattern)
            ]
            if not members:
                print(f"No {pattern} members in {archive}, keeping the archive", flush=True)
                continue
            print(f"Extracting {len(members)} file(s) from {archive}", flush=True)
            for member in members:
                # flatten: split_pepXML.py globs a single directory level
                target = os.path.join(os.path.dirname(archive), posixpath.basename(member))
                with z.open(member) as src, open(target, "wb") as dst:
                    shutil.copyfileobj(src, dst)
        os.remove(archive)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("credentials", help="yaml file with address/username/password/port")
parser.add_argument("remote", help="remote file, directory, or glob (in the last component)")
parser.add_argument("local", help="local file or directory to write")
parser.add_argument(
    "--unzip",
    metavar="GLOB",
    help="extract members matching GLOB out of any downloaded .zip archive and delete the archive",
)
args = parser.parse_args()

with open(args.credentials, "r") as inFile:
    cred = yaml.load(inFile, Loader=yaml.CLoader)

print(f"Connecting to {cred['address']}...", flush=True)
remote = Remote(cred)
location = remote.resolve(args.remote)

if remote.size(location) is None:
    print(f"Downloading directory {location} -> {args.local}", flush=True)
    downloaded = remote.download_tree(location, args.local)
else:
    os.makedirs(os.path.dirname(args.local) or ".", exist_ok=True)
    print(f"Downloading {location} -> {args.local}", flush=True)
    remote.download(location, args.local)
    downloaded = [args.local]

remote.ftp.quit()

if args.unzip:
    extract_zips(downloaded, args.unzip)
