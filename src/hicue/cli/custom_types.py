import click

from .imports import *

# accepts a list of cool files or a file containing a list of cool, one per line
class coolType(click.ParamType):
    name="cool"

    def is_cool(self, file):
        try:
            return cooler.fileops.is_cooler(file)
        except:
            return False

    def convert(self, value, param, ctx):
        if self.is_cool(value):
            return [value]
        if not os.path.isfile(value):
            files = value.split(',')
            for file in files:
                if not self.is_cool(file):
                    self.fail(f"{file} is not a valid cool file", param, ctx)
            return files
        else:
            files = []
            with open(value, 'r') as f:
                for line in f.readlines():
                    file = line.replace('\n', '').replace(' ', '')
                    if not self.is_cool(file):
                        self.fail(f"{file} is not a valid cool file", param, ctx)
                    else:
                        files.append(file)
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

    def convert(self, value, param, ctx):
        return value
    

COOL = coolType()
INT_LIST = IntListType()
STR_LIST = StrListType()
POSITION_FILE = PositionFileType()