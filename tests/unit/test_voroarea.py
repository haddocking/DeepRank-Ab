"""Unit tests for pure parsing helpers in src/tools/VoroArea.py.

Note: `VoronotaAreas.get_contact_areas` calls `self.get_atom_key`, which is not
defined anywhere in the class (no such method/staticmethod exists). This looks
like a real bug (dead code path / missing method) rather than something to
test around -- flagged for the user, not covered here.
"""

from src.tools.VoroArea import VoronotaAreas


# ---------------------------------------------------------------------------
# load_contacts
# ---------------------------------------------------------------------------

def test_load_contacts_parses_lines():
    raw = b"0 0 35.440\n0 1 15.908\n1 1 16.448\n"
    contacts = VoronotaAreas.load_contacts(raw)
    assert contacts == {0: {0: 35.440, 1: 15.908}, 1: {1: 16.448}}


def test_load_contacts_skips_blank_lines():
    raw = b"0 0 1.0\n\n\n1 1 2.0\n"
    contacts = VoronotaAreas.load_contacts(raw)
    assert contacts == {0: {0: 1.0}, 1: {1: 2.0}}


def test_load_contacts_empty_input():
    assert VoronotaAreas.load_contacts(b"") == {}


# ---------------------------------------------------------------------------
# load_balls
# ---------------------------------------------------------------------------

def test_load_balls_parses_lines():
    raw = (
        b"28.888 9.409 52.301 1.7 # 1 A 2 SER N\n"
        b"27.638 10.125 52.516 1.9 # 2 A 2 SER CA\n"
    )
    balls = VoronotaAreas.load_balls(raw)
    assert balls == {
        0: ("1", "A", "2", "SER", "N"),
        1: ("2", "A", "2", "SER", "CA"),
    }


def test_load_balls_skips_blank_lines():
    raw = b"28.888 9.409 52.301 1.7 # 1 A 2 SER N\n\n"
    balls = VoronotaAreas.load_balls(raw)
    assert len(balls) == 1
    assert balls[0] == ("1", "A", "2", "SER", "N")


def test_load_balls_empty_input():
    assert VoronotaAreas.load_balls(b"") == {}


# ---------------------------------------------------------------------------
# get_voronota_executable
# ---------------------------------------------------------------------------

def test_get_voronota_executable_path():
    exe = VoronotaAreas.get_voronota_executable()
    assert exe.endswith("src/tools/voronota/voronota")
