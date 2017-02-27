import subprocess

# Store version number in separate module to be accessible from all three: setup.py, gui.py, and pymfix.py
# version numbering from PEP 440
__version__ = "17.1dev"

def get_git_revision_short_hash():
    """Try to get the current git hash"""
    try:
        git_hash = subprocess.check_output(['git', 'describe', '--always']).strip().decode('utf-8')
    except:
        git_hash = None

    return git_hash
