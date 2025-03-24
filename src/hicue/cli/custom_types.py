import click

from .imports import *

class coolType(click.ParamType):
    name="cool"

    def convert(self, value, param, ctx):
        if cooler.fileops.is_cooler(value):
            return value
        else:
            self.fail(f"{value} is not a valid cool file", param, ctx)

class IntListType(click.ParamType):
    name="int_list"

    def convert(self, value, param, ctx):
        ints = value.split(',')
        for i in ints:
            if not isinstance(int(i), int):
                self.fail(f"{value} is not a comma-separated integer list.")
            if int(i) < 0:
                self.fail(f"All integers provided in {param.name} must be positive integers.")
        return value
    
class StrListType(click.ParamType):
    name="str_list"

    def convert(self, value, param, ctx):
        strs = value.split(',')
        for s in strs:
            if len(s) == 0:
                self.fail(f"Empty strings. Comma must separate two distinct values; if a single string is passed, no comma must be found in {param.name}.")
        return value
    
class PositionFileType(click.ParamType):
    name="position file"

    # TODO define PositionFileType
    def convert(self, value, param, ctx):
        return value
    

COOL = coolType()
INT_LIST = IntListType()
STR_LIST = StrListType()
POSITION_FILE = PositionFileType()