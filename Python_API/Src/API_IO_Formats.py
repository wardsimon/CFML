# **************************************************************************
#
# CrysFML API IO Format
#
# @file      Src/API_IO_Formats.py
# @brief     IO Formats based on CrysFML
#
# @homepage  https://code.ill.fr/scientific-software/crysfml
# @license   GNU LGPL (see LICENSE)
# @copyright Institut Laue Langevin 2020-now
# @authors   Scientific Computing Group at ILL (see AUTHORS)
#
# **************************************************************************

import os

import CFML_api.crysfml_api
import CFML_api.API_Crystal_Metrics
import CFML_api.API_Crystallographic_Symmetry
import CFML_api.API_Atom_TypeDef
import CFML_api.FortranBindedClass

class JobInfo(CFML_api.FortranBindedClass):
    """ A class used to handle the job informations, read from a CFL File

    ...
    Attributes
    ----------
    string_array:

    Methods
    -------

    """
    def __init__(self, string_array=None):
        CFML_api.FortranBindedClass.__init__(self)
        if string_array:
            dict = CFML_api.crysfml_api.IO_Formats_jobinfo_from_CIF_string_array(string_array)
            self._set_fortran_address(dict["JobInfo"])
            

    def __del__(self):
        address = self.get_fortran_address()
        if address:
            CFML_api.crysfml_api.IO_Formats_del_jobinfo(self.get_fortran_address())
    
    @property
    def title(self):
        """
        Identification of the job
        """
        return CFML_api.crysfml_api.IO_Formats_get_title(self.get_fortran_address())["title"]
    
    @property
    def num_phases(self):
        """
        Number of phases
        """
        return CFML_api.crysfml_api.IO_Formats_get_num_phases(self.get_fortran_address())["num_phases"]
    
    @property
    def num_patterns(self):
        """
        Number of patterns
        """
        return CFML_api.crysfml_api.IO_Formats_get_num_patterns(self.get_fortran_address())["num_patterns"]
    
    @property
    def num_cmd(self):
        """
        Number of command lines
        """
        return CFML_api.crysfml_api.IO_Formats_get_num_cmd(self.get_fortran_address())["num_cmd"]
    
    @property
    def pattern_types(self, indx=None):
        """
        String with the pattern type of pattern indx in [0, num_patterns-1]
        """
        if indx:
            key=indx
        else:
            key=0
        return CFML_api.crysfml_api.IO_Formats_get_patt_typ(self.get_fortran_address(), key+1)["patt_typ"]
    
    @property
    def phase_names(self, indx=None):
        """
        String with the name of the phase indx in [0, num_phases-1]
        """
        if indx:
            key=indx
        else:
            key=0
        return CFML_api.crysfml_api.IO_Formats_get_phas_nam(self.get_fortran_address(), key+1)["phas_nam"]
    
    @property
    def cmd(self, indx=None):
        """
        Command lines: text for actions in [0, num_cmd-1]
        """
        if indx:
            key=indx
        else:
            key=0
        return CFML_api.crysfml_api.IO_Formats_get_cmd(self.get_fortran_address(),key+1)["cmd"]
    
    @property
    def range_stl(self,indx=None):
        """
        Range in sinTheta/Lambda of phase indx
        """
        if indx:
            key=indx
        else:
            key=0
        mina = CFML_api.crysfml_api.IO_Formats_get_range_stl(self.get_fortran_address(),key+1)["min"]
        maxb = CFML_api.crysfml_api.IO_Formats_get_range_stl(self.get_fortran_address(),key+1)["max"]
            
        return mina, maxb
    
    @property
    def range_q(self,indx=None):
        """
        Range in 4pi*sinTheta/Lambda of phase indx
        """
        if indx:
            key=indx
        else:
            key=0
        mina = CFML_api.crysfml_api.IO_Formats_get_range_q(self.get_fortran_address(),key+1)["min"]
        maxb = CFML_api.crysfml_api.IO_Formats_get_range_q(self.get_fortran_address(),key+1)["max"]
            
        return mina, maxb
        
    
    @property
    def range_d(self, indx=None):
        """
        Range in d-spacing of phase indx
        """
        if indx:
            key=indx
        else:
            key=0
        mina = CFML_api.crysfml_api.IO_Formats_get_range_d(self.get_fortran_address(),key+1)["min"]
        maxb = CFML_api.crysfml_api.IO_Formats_get_range_d(self.get_fortran_address(),key+1)["max"]
            
        return mina, maxb
        
    
    @property
    def range_2theta(self, indx=None):
        """
        Range in 2theta-spacing of phase indx
        """
        if indx:
            key=indx
        else:
            key=0
        mina = CFML_api.crysfml_api.IO_Formats_get_range_2theta(self.get_fortran_address(),key+1)["min"]
        maxb = CFML_api.crysfml_api.IO_Formats_get_range_2theta(self.get_fortran_address(),key+1)["max"]
            
        return mina, maxb
        
    
    @property
    def range_energy(self, indx=None):
        """
        Range in Energy of phase indx
        """
        if indx:
            key=indx
        else:
            key=0
        mina = CFML_api.crysfml_api.IO_Formats_get_range_energy(self.get_fortran_address(),key+1)["min"]
        maxb = CFML_api.crysfml_api.IO_Formats_get_range_energy(self.get_fortran_address(),key+1)["max"]
            
        return mina, maxb
        
    
    @property
    def lambdas(self, indx=None):
        """
        Lambda1, Lambda2 of phase indx
        """
        if indx:
            key=indx
        else:
            key=0
        lambda1 = CFML_api.crysfml_api.IO_Formats_get_lambda(self.get_fortran_address(),key+1)["lambda1"]
        lambda2 = CFML_api.crysfml_api.IO_Formats_get_lambda(self.get_fortran_address(),key+1)["lambda2"]
            
        return lambda1, lambda2
 
    @property
    def lambda_ratio(self):
        """
        ratio lambda2/lambda1 
        """
        return CFML_api.crysfml_api.IO_Formats_get_ratio(self.get_fortran_address())["ratio"]
    
    @property
    def d_to_tof_1(self):
        """
        d-to-TOF coefficient
        """
        return CFML_api.crysfml_api.IO_Formats_get_dtt1(self.get_fortran_address())["dtt1"]
    
    @property
    def d_to_tof_2(self):
        """
        d-to-TOF coefficient
        """
        return CFML_api.crysfml_api.IO_Formats_get_dtt2(self.get_fortran_address())["dtt2"]
    

class CIFFile():
    """ A class used to handle CIF Files

    ...
    Attributes
    ----------
    filename : string

    Methods
    -------
    get_as_string_array

    """
    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str
           Name of the CIF File
        """
        self.__load_cif_file(filename)
    
    @property
    def filename(self):
        """
        Name of the CIF File (str)
        """
        return self.__filename
    
    @property
    def cell(self):
        """
        Cell (CFML_api.Cell type)
        """
        return self.__cell
    
    @property
    def space_group(self):
        """
        Space Group (CFML_api.SpaceGroup type)
        """
        return self.__space_group
    
    @property
    def atom_list(self):
        """
        List of atoms (CFML_api.AtomList type)
        """
        return self.__atom_list

    @property
    def job_info(self):
        """
        Info on the type of job (CFML_api.Job_info type)
        """
        return self.__job_info
    
    @filename.setter
    def filename(self, filename):
        """
        Setter for the name of the cif file
        The setter reloads the cell, space group and atom list
        """
        self.__load_cif_file(filename)    

    def get_as_string_array(self):
        """
        Returns the CIF File as an array of string
        """
        ret = []
        with open(self.__filename, 'r') as file:
            for line in file.readlines():
                ret.append(line[:-1])
        return ret

    def __load_cif_file(self, filename):
        """
        Internal method
        Loads the CIF File
        """
        self.__filename = filename

        # Check if file exists
        file_exists = False
        if filename[-4:] == ".cif" or filename[-4:] == ".cfl":
            if os.path.exists(filename):
                file_exists = True
                if filename[-4:] == ".cif":
                    mode = "CIF"
                else:
                    mode = "CFL"
        elif os.path.exists(filename+".cfl"):
                filename = filename + ".cfl"
                file_exists = True
                mode = "CFL"
        elif os.path.exists(filename+".cif"):
                filename = filename + ".cif"
                file_exists = True
                mode = "CIF"

        if not file_exists:
            raise IOError("No file with name: " + filename)
        else:
                
            dict = CFML_api.crysfml_api.IO_Formats_readn_set_xtal_structure(self.__filename, mode)
            print('mode',mode)
            self.__cell = CFML_api.API_Crystal_Metrics.Cell.from_fortran_address(dict["Cell"])
            self.__space_group = CFML_api.API_Crystallographic_Symmetry.SpaceGroup.from_fortran_address(dict["SpG"])
            self.__atom_list = CFML_api.API_Atom_TypeDef.AtomList.from_fortran_address(dict["A"])
            self.__job_info = JobInfo.from_fortran_address(dict["JobInfo"])
            #self.__job_info = CFML_api.API_IO_Formats.JobInfo.from_fortran_address(dict["JobInfo"])
