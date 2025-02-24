import subprocess
import platform
import stat
import sys
import os

def detect_executable():
    system = platform.system()
    machine = platform.machine().lower()

    # select correct muscle file based on OS and architecture
    if system == "Windows":
        return "executables/muscle-win64.v5.3.exe"
    elif system == "Linux":
        if "aarch64" in machine or "arm64" in machine:
            os.chmod("executables/muscle-aarch64.v5.3", stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
            return "executables/muscle-aarch64.v5.3"
        elif "x86_64" in machine or "amd64" in machine:
            os.chmod("executables/muscle-linux-x86.v5.3", stat.S_IRWXU | stat.S_IRWXG | stat.S_IRWXO)
            return "executables/muscle-linux-x86.v5.3"
    elif system == "Darwin":  # macOS
        if "arm64" in machine:
            return "executables/muscle-osx-arm64.v5.3"
        elif "x86_64" in machine:
            return "executables/muscle-osx-x86.v5.3"
    else:
        print(f"Unsupported system: {system} with architecture: {machine}")
        sys.exit(1)

def run_muscle(parameters):
    executable = detect_executable()

    # make sure muscle file exists
    if not os.path.isfile(executable):
        print(f"Executable not found: {executable}")
        sys.exit(1)

    # for macOS/Linux, ensure the file has execution permission
    # if platform.system() != "Windows":
    #     os.chmod(executable, 0o755)

    subprocess.run([executable] + parameters)

if __name__ == "__main__":

    in_file = sys.argv[1]
    out_file = sys.argv[2]

    os.makedirs(os.path.dirname(out_file), exist_ok=True)

    parameters = ["-align", in_file, "-output", out_file]
    run_muscle(parameters)

