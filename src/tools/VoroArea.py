import subprocess
from pathlib import Path
import numpy as np

class VoronotaAreas():
    def __init__(self, pdb_path):
        self.voronota_exec = self.get_voronota_executable()
        self.contact_areas = self._get_voronota_contacts(pdb_path)

    def _get_voronota_contacts(self, pdb_path) -> dict[tuple, dict[tuple, float]]:
        balls_outputs, contacts_outputs = self.run_voro_contacts(pdb_path, self.voronota_exec)
        # Load data
        balls_data = self.load_balls(balls_outputs)
        contacts_data = self.load_contacts(contacts_outputs)
        # map both together into a single dict
        contact_areas = {
            balls_data[at1_index]: {
                balls_data[at2_index]: contact_area
                for at2_index, contact_area in at1_contact_areas.items()
                }
            for at1_index, at1_contact_areas in contacts_data.items()
        }
        return contact_areas

    @staticmethod
    def load_contacts(_contacts_outputs) -> dict[int, dict[int, float]]:
        """
    In "contacts.txt" file the line format is "b1 b2 area".
    The first two numbers (b1 and b2) are numbers of atomic records in "balls.txt", starting from 0.
    If b1 does not equal b2, then the 'area' value is the area of contact between atoms b1 and b2.
    If b1 equals b2, then the 'area' value is the solvent-accessible area of atom b1.
    For example, below is a part of some possible "contacts.txt":

    0 0 35.440
    0 1 15.908
    0 2 0.167
    0 3 7.025
    0 4 7.021
    0 5 0.624
    0 23 2.849
    0 25 0.008
    0 26 11.323
    0 1454 0.021
    1 1 16.448
    1 2 11.608
        """
        lines = _contacts_outputs.decode("utf-8")
        contacts_data = {}
        for index, _ in enumerate(lines.split("\n")):
            striped_ = _.strip()
            if striped_ == "":
                continue
            s_ = striped_.split(" ")
            at1_index = int(s_[0])
            at2_index = int(s_[1])
            area = float(s_[2])
            at1_areas = contacts_data.setdefault(at1_index, {})
            at1_areas[at2_index] = area
        return contacts_data

    @staticmethod
    def load_balls(_balls_data) -> dict[int, tuple]:
        """
    In "balls.txt" the line format is "x y z r # comments".
    The first four values (x, y, z, r) are atomic ball coordinates and radius.
    Comments are not needed for further calculations, they are to assist human readers.
    For example, below is a part of some possible "balls.txt":

    28.888 9.409 52.301 1.7 # 1 A 2 SER N
    27.638 10.125 52.516 1.9 # 2 A 2 SER CA
    26.499 9.639 51.644 1.75 # 3 A 2 SER C
    26.606 8.656 50.915 1.49 # 4 A 2 SER O
    27.783 11.635 52.378 1.91 # 5 A 2 SER CB
    27.69 12.033 51.012 1.54 # 6 A 2 SER OG
        """
        lines = _balls_data.decode("utf-8")
        balls_data = {}
        for index, _ in enumerate(lines.split("\n")):
            if _.strip() == "":
                continue
            atom_dt = _.split("#")[-1].strip().split(" ")
            desired_at_dt = tuple(atom_dt[:5])
            balls_data[index] = desired_at_dt
        return balls_data

    @staticmethod
    def get_voronota_executable() -> str:
        # Absolute path to the directory of this file
        script_dir = Path(__file__).resolve().parent

        # DeepRank-Ab root directory
        root_dir = script_dir.parent

        # Tools directories
        voronota_exec = root_dir / "tools" / "voronota" / "voronota"

        return str(voronota_exec)

    # TODO: function too complex, refactor
    @staticmethod
    def run_voro_contacts(pdb_fpath, voronota_exec) -> tuple[bytes, bytes]:
        # Read input file as bytes
        with open(pdb_fpath, "rb") as fin:
            pdb_lines = fin.read()

        # Setup the two subproces commands
        balls = subprocess.Popen(
            [voronota_exec, "get-balls-from-atoms-file"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            )
        contacts = subprocess.Popen(
            [voronota_exec, "calculate-contacts"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            )

        # Run commands with inputs
        #./voronota get-balls-from-atoms-file < input.pdb > balls.txt
        balls_outputs = balls.communicate(input=pdb_lines)[0]
        #./voronota calculate-contacts < balls.txt > contacts.txt
        contacts_outputs = contacts.communicate(input=balls_outputs)[0]

        # Return outputs of both commands
        return balls_outputs, contacts_outputs

    def get_contact_areas(self, atoms1, atoms2) -> np.float64:
        areas = []
        for at1 in atoms1:
            at1k = self.get_atom_key(at1)
            for at2 in atoms2:
                # Set the contact between two atoms from same residue to 0
                if at1.residue == at2.residue:
                    areas.append(0)
                    continue

                # Get key for second atom
                at2k = self.get_atom_key(at2)
                try:
                    area = self.contact_areas[at1k][at2k]
                except KeyError:
                    try:
                        area = self.contact_areas[at2k][at1k]
                    except KeyError:
                        area = 0
                areas.append(area)
        return np.sum(areas)
