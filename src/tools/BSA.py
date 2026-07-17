import os
import numpy as np
from pdb2sql.interface import interface

try:
    import freesasa

except ImportError:
    print('Freesasa not found')


# TODO: function too complex, refactor
class BSA(object):

    def __init__(self, pdb_data, sqldb=None, chainA='A', chainB='B'):
        '''Compute the burried surface area feature

        Freesasa is required for this feature.

        https://freesasa.github.io

        >>> wget http://freesasa.github.io/freesasa-2.0.2.tar.gz
        >>> tar -xvf freesasa-2.0.3.tar.gz
        >>> cd freesasa
        >>> ./configure CFLAGS=-fPIC (--prefix /home/<user>/)
        >>> make
        >>> make install

        Since release 2.0.3 the python bindings are separate module
        >>> pip install freesasa

        Args :
            pdb_data (list(byte) or str): pdb data or filename of the pdb
            sqldb (pdb2sql.interface instance or None, optional) if the sqldb is None the sqldb will be created
            chainA (str, optional): name of the first chain
            chainB (str, optional): name of the second chain

        Example :

        >>> bsa = BSA('1AK4.pdb')
        >>> bsa.get_structure()
        >>> bsa.get_contact_residue_sasa()
        >>> bsa.sql.close()

        '''

        self.pdb_data = pdb_data
        if sqldb is None:
            self.sql = interface(pdb_data)
        else:
            self.sql = sqldb
        self.chains_label = [chainA, chainB]

        freesasa.setVerbosity(freesasa.nowarnings)

    def get_structure(self):
        """This method prepares three Freesasa “Structure” objects:
        The full complex

        Chain A alone

        Chain B alone

    """
        """Get the pdb structure of the molecule."""
# 1. FULL COMPLEX
        # we can have a str or a list of bytes as input
        if isinstance(self.pdb_data, str):
            # If you passed a filename, it just hands that straight to freesasa.Structure(pdbfile), which parses it itself.
            self.complex = freesasa.Structure(self.pdb_data)
        else:
            # if you passed a list of bytes, it will parse it itself. 
            # create empty structure
            self.complex = freesasa.Structure()
            # take the columns from the pdb data
            atomdata = self.sql.get(
                'name,resName,resSeq,chainID,x,y,z')
            # reformat and add the atoms to the structure
            for atomName, residueName, residueNumber, chainLabel, x, y, z in atomdata:
                atomName = '{:>2}'.format(atomName[0])
                self.complex.addAtom(
                    atomName, residueName, residueNumber, chainLabel, x, y, z)
        #############   COMPUTE THE BSA FOR THE COMPLEX   #############
        self.result_complex = freesasa.calc(self.complex)



# 2. CHAINS
        # we will store here the fresasa structure objects for each chain
        self.chains = {}
        # we will store here the result of the freesasa calc for each chain
        self.result_chains = {}
        # for each chain
        for label in self.chains_label:
            # empty structure
            self.chains[label] = freesasa.Structure()
            # adding the atoms to the structure
            atomdata = self.sql.get(
                'name,resName,resSeq,chainID,x,y,z', chainID=label)
            for atomName, residueName, residueNumber, chainLabel, x, y, z in atomdata:
                atomName = '{:>2}'.format(atomName[0])
                self.chains[label].addAtom(
                    atomName, residueName, residueNumber, chainLabel, x, y, z)
            self.result_chains[label] = freesasa.calc(
                self.chains[label])

    def get_contact_residue_sasa(self, cutoff=8.5):
        """Compute the feature value."""

        # value
        self.bsa_data = {}
        # where it is 
        self.bsa_data_xyz = {}

        # find contact residues
        res = self.sql.get_contact_residues(cutoff=cutoff)
        # save theorder of the chains
        keys = list(res.keys())
        # concatenates the two lists containing the residues of the interface in both chains
        res = res[keys[0]]+res[keys[1]]

        # for each residue in the interface
        # ex. r = ('A', 125)
        for r in res:

            # sselect the residue and the chain
            select_str = ('res, (resi %d) and (chain %s)' %
                          (r[1], r[0]),)
            # free sasa tells us the saea  (strs)
            asa_complex = freesasa.selectArea(
                select_str, self.complex, self.result_complex)['res']

            # same when the residue is isolated
            select_str = ('res, resi %d' % r[1],)
            asa_unbound = freesasa.selectArea(
                select_str, self.chains[r[0]], self.result_chains[r[0]])['res']

            # define the bsa (area unbound - area complex)
            bsa = asa_unbound-asa_complex

            # define the xyz key : (chain,x,y,z)
            chain = {'A': 0, 'B': 1}[r[0]] # assign 0 to A and 1 to B
            # get the average xyz of the residue to get the
            # geometric center of the residue
            xyz = np.mean(self.sql.get(
                'x,y,z', resSeq=r[1], chainID=r[0]), 0)
            
            # [0, 12.3, -7.8, 5.6]
            xyzkey = tuple([chain]+xyz.tolist())

            # put the data in dict
            self.bsa_data[r] = [bsa]
            self.bsa_data_xyz[r] = xyzkey

