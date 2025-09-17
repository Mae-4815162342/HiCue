import click

from .imports import *

from hicue.parser import *

# accepts a list of cool files or a file containing a list of cool, one per line
class coolType(click.ParamType):
    name="cool"

    def is_cool(self, file):
        try:
            return cooler.fileops.is_cooler(file)
        except:
            return False
        
    def is_mcool(self, file):
        try:
            return cooler.fileops.is_multires_file(file)
        except:
            return False

    def convert(self, value, param, ctx):
        if self.is_cool(value) or self.is_mcool(value):
            return [value]
        if not os.path.isfile(value):
            files = value.split(',')
            for file in files:
                if not self.is_cool(file) and not self.is_mcool(file):
                    self.fail(f"{file} is not a valid cool file", param, ctx)
            return files
        else:
            files = []
            with open(value, 'r') as f:
                for line in f.readlines():
                    file = line.replace('\n', '').replace(' ', '')
                    if not self.is_cool(file) and not self.is_mcool(file):
                        self.fail(f"{file} is not a valid cool file", param, ctx)
                    else:
                        files.append(file)
            return files

# accepts a pair of comma-separated cool files or a file containing a list of pairs of cool, one comma-separated pair per line
class coolPairType(click.ParamType):
    name="cool_pairs"

    def is_cool(self, file):
        try:
            return cooler.fileops.is_cooler(file)
        except:
            return False
    
    def is_cool_pair(self, value):
        splited_values = value.split(',')
        if len(splited_values) != 2:
            return False
        if not self.is_cool(splited_values[0]) or not self.is_cool(splited_values[1]):
            return False
        return True
        
    def is_same_bin_pair(self, value):
        splited_values = value.split(',')
        # checking the pair has the same binning
        cool1 = cooler.Cooler(splited_values[0])
        cool2 = cooler.Cooler(splited_values[1])
        if cool1.binsize != cool2.binsize:
            return False
        return True

    def convert(self, value, param, ctx):
        if self.is_cool_pair(value):
            if not self.is_same_bin_pair(value):
                self.fail(f"The files {value} do not have the same binning value.", param, ctx)
            return [value.split(',')]
        elif not os.path.isfile(value):
            self.fail(f"{value} is not a valid comma-separated pair of cool file, neither a file containing pairs of cool files.", param, ctx)
        else:
            files = []
            with open(value, 'r') as f:
                for line in f.readlines():
                    file_list = line.replace('\n', '').replace(' ', '')
                    if not self.is_cool_pair(file_list):
                        self.fail(f"{file_list} is not a valid comma-separated pair of cool file.", param, ctx)
                    if not self.is_same_bin_pair(file_list):
                        self.fail(f"The files {file_list} do not have the same binning value.", param, ctx)
                    else:
                        files.append(file_list.split(","))
            return files

class IntListType(click.ParamType):
    name="int_list"

    def convert(self, value, param, ctx):
        list = value.split(',')
        ints = []
        for i in list:
            if not isinstance(int(i), int):
                self.fail(f"{value} is not a comma-separated integer list.")
            if int(i) < 0:
                self.fail(f"All integers provided in {param.name} must be positive integers.")
            ints.append(int(i))
        return ints
    
class StrListType(click.ParamType):
    name="str_list"

    def convert(self, value, param, ctx):
        strs = value.split(',')
        for s in strs:
            if len(s) == 0:
                self.fail(f"Empty strings. Comma must separate two distinct values; if a single string is passed, no comma must be found in {param.name}.")
        return strs

class PositionFileType(click.ParamType):
    name="position file"

    #  accepts bed and gff files. Returns a tuple indicating the type and path of the file
    def convert(self, value, param, ctx):
        """Accepts gff, bed and bed2d files."""
        if not os.path.isfile(value):
            self.fail(f"{value} is not an existing file.")
        extension = value.split('.')[-1].lower()
        if extension not in ["gff", "bed", "bed2d"]:
            self.fail(f"{value} extension is not recognized (gff, bed, or bed2d expected)")
        return extension, value
    
class Position2dFileType(click.ParamType):
    name="position file"

    #  accepts bed and gff files. Returns a tuple indicating the type and path of the file
    def convert(self, value, param, ctx):
        if not os.path.isfile(value):
            self.fail(f"{value} is not an existing file.")
        try:
            parse_bed2d_file(value)
            return value
        except Exception as e:
            self.fail(f"{value} is not an valid 2d position file (bed2d format expected) : {e}")

class GffFileType(click.ParamType):
    name="gff file"

    #  accepts bed and gff files. Returns a tuple indicating the type and path of the file
    def convert(self, value, param, ctx):
        if not os.path.isfile(value):
            self.fail(f"{value} is not an existing file.")
        try:
            in_handle = open(value)
            recs = GFF.parse(in_handle)
            if len(list(recs)) == 0:
                raise Exception
            return value
        except Exception as e:
            self.fail(f"{value} is not an valid gff file : {e}")

class TrackFileType(click.ParamType):
    name="track file"

    #  accepts bw files
    def convert(self, value, param, ctx):
        try:
            bw = pyBigWig.open(value)
            if bw.isBigWig():
                return value
            else:
                self.fail(f"{value} provided for {param.name} is not in the bigWig format.")
        except:
            self.fail(f"{value} provided for {param.name} is not in the bigWig format.")



COOL = coolType()
COOL_PAIR = coolPairType()
INT_LIST = IntListType()
STR_LIST = StrListType()
POSITION_FILE = PositionFileType()
POSITION2D_FILE = Position2dFileType()
TRACK_FILE = TrackFileType()
GFF_FILE = GffFileType()