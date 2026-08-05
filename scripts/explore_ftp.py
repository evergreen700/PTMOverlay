import ftplib
import yaml
import sys

credentials = sys.argv[1]
path = sys.argv[2] if len(sys.argv) > 2 else "."

with open(credentials) as f:
    cred = yaml.load(f, Loader=yaml.CLoader)

ftp = ftplib.FTP_TLS()
ftp.connect(cred["address"], cred.get("port", 21))
ftp.login(user=cred.get("username", "anonymous"), passwd=cred.get("password", ""))

print(f"PWD: {ftp.pwd()}")
try:
    ftp.size(path)
    lines = []
    ftp.retrlines(f"RETR {path}", lines.append)
    print("\n".join(lines))
except ftplib.error_perm:
    ftp.dir(path)
ftp.quit()
