#!/usr/bin/env python
import sys, re, subprocess, pytest
from tempfile import TemporaryDirectory


# Change into directory where test_recipy.py lives
# As suggested by https://stackoverflow.com/questions/62044541/change-pytest-working-directory-to-test-case-directory
@pytest.fixture(autouse=True)
def change_test_dir(request, monkeypatch):
    monkeypatch.chdir(request.fspath.dirname)

def run(command):
    """Runs command, returns tuple of exit code, output"""
    print(f"running: {command}")
    try:
        return 0, subprocess.check_output(command, shell=True).strip().decode()
    except subprocess.CalledProcessError as exc:
        return exc.returncode, exc.output.strip().decode()

def verify_output(output, *regexes):
    """Returns true if the output contains lines matching the regexes in sequence (possibly with other lines in between)"""
    regexes = list(regexes[::-1])
    for line in output.split("\n"):
        if regexes and re.search(regexes[-1], line):
            regexes.pop()
    if regexes:
        print("Error, the following regexes did not match the output:")
        for regex in regexes:
            print(f"  {regex}")
        return False
    return True

def test_convert(tmpdir: str = None):
    if tmpdir is None:
        tmpdir_obj = TemporaryDirectory()
        tmpdir = tmpdir_obj.name
        print(f"temporary output dir is {tmpdir}")

    print("Test conversion")

    retcode, _ = run(f"""tigger-convert 3C147-HI6.refmodel.lsm.html {tmpdir}/output.txt -f """
        """ --output-format "name ra_d dec_d i q u v i q u v spi rm emaj_s emin_s pa_d freq0" """)
    assert retcode == 0

    retcode, _ = run(f"""tigger-convert 3C147-HI6.refmodel.lsm.html {tmpdir}/output.recentred.txt -f """
        """--center 85.5deg,49.9deg --rename --output-format "name ra_d dec_d i q u v i q u v spi rm emaj_s emin_s pa_d freq0" """)
    assert retcode == 0

    print("Checking reference LSM")
    retcode, _ = run(f"""diff 3C147-HI6.refmodel.reference.txt {tmpdir}/output.txt""")
    assert retcode == 0

    print("Checking recentred reference LSM")

    retcode, _ = run(f"""diff 3C147-HI6.refmodel.recentred.reference.txt {tmpdir}/output.recentred.txt""")
    assert retcode == 0

    print("Test reverse conversion to .lsm.html models")
    retcode, _ = run(f"""tigger-convert {tmpdir}/output.txt {tmpdir}/output.lsm.html -f""")
    assert retcode == 0

    print("Test .gaul conversions")
    retcode, _ = run(f"""tigger-convert deep4.gaul {tmpdir}/deep4.lsm.html -f""")
    assert retcode == 0
    retcode, _ = run(f"""tigger-convert deep4.gaul {tmpdir}/deep4.txt -f --output-format "name ra_d dec_d i q u v spi rm emaj_s emin_s pa_d freq0" """)
    assert retcode == 0
    retcode, _ = run(f"""diff deep4.reference.txt {tmpdir}/deep4.txt""")
    assert retcode == 0
    
    print("Test .AIPSCC conversions")
    retcode, _ = run(f"""gunzip <3C147-L-A-CLEAN.fits.gz >{tmpdir}/3C147-L-A-CLEAN.fits""")
    assert retcode == 0
    retcode, _ = run(f"""tigger-convert {tmpdir}/3C147-L-A-CLEAN.fits {tmpdir}/3C147-L-A-CLEAN.fits.lsm.html -f""")
    assert retcode == 0
    retcode, _ = run(f"""tigger-convert {tmpdir}/3C147-L-A-CLEAN.fits.lsm.html {tmpdir}/3C147-L-A-CLEAN.txt -f --output-format "name ra_d dec_d i q u v" """)
    assert retcode == 0
    retcode, _ = run(f"""zdiff 3C147-L-A-CLEAN.txt.gz {tmpdir}/3C147-L-A-CLEAN.txt""")
    assert retcode == 0

    print("Testing tigger-restore and tigger-make-brick")
    retcode, _ = run(f"""cp 3C147tmp.fits {tmpdir}""")
    assert retcode == 0
    retcode, _ = run(f"""tigger-make-brick 3C147-HI6.refmodel.lsm.html {tmpdir}/3C147tmp.fits """)
    assert retcode == 0
    retcode, _ = run(f"""tigger-restore -f 3C147tmp.fits 3C147-HI6.refmodel.lsm.html {tmpdir}/restored.fits""")
    assert retcode == 0
    retcode, _ = run(f"""tigger-make-brick 3C147-HI6.refmodel.lsm.html {tmpdir}/3C147tmp.fits {tmpdir}/brickmodel.lsm.html -f""")
    assert retcode == 0
    retcode, _ = run(f"""tigger-restore -f 3C147tmp.fits 3C147-HI6.refmodel.lsm.html {tmpdir}/restored1.fits""")
    assert retcode == 0
    retcode, _ = run(f"""zdiff {tmpdir}/restored1.fits restored1-reference.fits.gz""")
    assert retcode == 0

    print("Testing tigger-convert add-brick")
    retcode, _ = run(f"""tigger-convert 3C147-HI6.refmodel.lsm.html {tmpdir}/add-brick.lsm.html -f """
        """ --add-brick BRICK:3C147tmp-zoom.fits """)
    assert retcode == 0

    print("Testing tigger-restore")
    retcode, _ = run(f"""tigger-restore -f 3C147tmp.fits {tmpdir}/add-brick.lsm.html {tmpdir}/restored2.fits""")
    assert retcode == 0
    retcode, _ = run(f"""zdiff {tmpdir}/restored2.fits restored2-reference.fits.gz""")
    assert retcode == 0

    print("Test tigger-tag")
    retcode, _ = run(f"""tigger-tag 3C147-HI6.refmodel.lsm.html 'r<0.5d' inner=1 -o {tmpdir}/tmp.lsm.html -f """)
    assert retcode == 0

    print("Test .lsm to ds9 region file conversion")
    retcode, _ = run(f"""tigger-convert 3C147-HI6.refmodel.lsm.html {tmpdir}/3C147-HI6.refmodel.lsm.reg -f """)
    assert retcode == 0

    return tmpdir

if __name__ == '__main__':
    if len(sys.argv) > 1:
        tmpdir = sys.argv[1]
    else:
        tmpdir = None

    tmpdir = test_convert(tmpdir)

    print(f"Output directory was {tmpdir}")