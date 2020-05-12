import re

class Config_file_parser():
    """ This class is used for parsing the configuration file "conf.txt".  """
    def __init__(self, filepath=None, verbose=1):
        self.filepath = filepath
        self.regex_dict = dict()
        self.key_type = dict()
        self.key_ownprocess = dict()
        self.data = self.init_data()
        self.verbose = verbose

    def init_data(self):
        """ Initialize data with default values for example """
        return {"todo": None, "unit_cells": 1, "plotphase":False}

    def add_pattern(self, key=None, key_type=str, regex='default', ownprocess=None, **kwargs):
        """
        key:        string identifying the variable to parse in the config file
        key_type:   the type which the value behind key (initially text) must be
                    converted to.
        regex:      the regular expression for finding key and its value when 
                    parsing the config file.
        ownprocess:
        """
        if 'regex_dict' in kwargs:
            if isinstance(kwargs['regex_dict'], dict):
                re_dict = dict()
                type_dict = dict()
                for k,it in kwargs['regex_dict'].items():
                    re_dict[k] = it[0]
                    type_dict[k] = it[1]
                self.regex_dict.update(re_dict)
                self.key_type.update(type_dict)
        if key != None:
            if regex == 'default': regex = key + '=(?P<' + key + '>.*)'
            self.regex_dict[key] = re.compile(regex)
            self.key_type[key] = key_type
            if ownprocess != None:
                self.key_ownprocess[key] = ownprocess

    def _parse_line(self, line):
        """ Read a line and parse it with the regex dict.  """
        for key, rx in self.regex_dict.items():
            match = rx.match(line)
            if match:
                return key, match
        return None, None

    def _process_line(self, line):
        """ process a line with the regex dict.  """
        key, match_object = self._parse_line(line)
        if key == None or match_object == None: return

        if self.key_type[key] == float:
            self.data[key] = float(match_object.group(key))
        elif self.key_type[key] == int:
            self.data[key] = int(match_object.group(key))
        elif self.key_type[key] == str:
            self.data[key] = match_object.group(key)
        elif self.key_type[key] == bool:
            if match_object.group(key).lower() == 'true': self.data[key] = True
            else: self.data[key] = False
        elif self.key_type[key] == None:
            if key in self.key_ownprocess:
                self.key_ownprocess[key](key, match_object, self.data)
        else:
            print('WARNING: the type ' + self.key_type[key] + ' is not yet treated')

    def process_file(self):
        """ process each line """
        with open(self.filepath, 'r') as file_object:
            line = file_object.readline()
            while line:
                if self.verbose == 1:
                    print('Processing:' + str(line.encode('unicode-escape')))
                line = line.replace(' ','')
                if line[0] != '#': # do not process commentary
                    self._process_line(line)
                line = file_object.readline()
        return self.data





def retrieve_layer(key, match_object, data):
    """ Special regex for layers"""
    if key == 'layers':
        if key not in data:
            data[key] = {}
        idname = match_object.group('layer_id')
        name = match_object.group('layer_name')
        thickness = float(match_object.group('layer_thickness'))
        if name=='meta':
            Mat = match_object.group('Mat')
            Inc = match_object.group('Inc')
            rmean = float(match_object.group('rmean'))
            phi = float(match_object.group('phi'))
            poly = float(match_object.group('poly'))
            data[key][idname] = (name, thickness, Mat, Inc, rmean, phi, poly)
        else:
            data[key][idname] = (name, thickness)





def read_config_file(filepath, v):
    #interpreting config file with the object Config_file_parser
    config_reader = Config_file_parser(filepath, verbose=v)
    #set up the protocol to retrieve the parameters from the config file
    config_reader.add_pattern('f_min', float, regex='default')
    config_reader.add_pattern('f_max', float, regex='default')
    config_reader.add_pattern('f_fix', float, regex='default')
    config_reader.add_pattern('f_num', int, regex='default')
    config_reader.add_pattern('theta_min', float, regex='default')
    config_reader.add_pattern('theta_max', float, regex='default')
    config_reader.add_pattern('theta_fix', float, regex='default')
    config_reader.add_pattern('theta_num', int, regex='default')
    # config_reader.add_pattern('f_loop', bool, regex='default')
    # config_reader.add_pattern('nb_layers', int, regex='default')
    config_reader.add_pattern('unit_cells', int, regex='default')
    config_reader.add_pattern('todo', str, regex='default')
    # regex_layers = r'layer\#(?P<layer_id>\d+)={(?P<layer_name>.*):(?P<layer_thickness>.*)}'
    # regex_layers = r'layer\#(?P<layer_id>\d+)\s?=\s?{(?P<layer_name>[^ ,:]*)\s?:\s?(?P<layer_thickness>\d*)(\s?,\s?Mat\s?:\s?(?P<Mat>.*),\s?Inc\s?:\s?(?P<Inc>.*)\s?,\s?rmean\s?:\s?(?P<rmean>.*)\s?,\s?phi\s?:\s?(?P<phi>.*)\s?,\s?poly\s?:\s?(?P<poly>.*).*}|})'
    regex_layers = r'layer\#(?P<layer_id>\d+)\s?=\s?{(?P<layer_name>[^ ,:]*)\s?:\s?(?P<layer_thickness>[-+]?[0-9]*\.?[0-9]*)(\s?,\s?Mat\s?:\s?(?P<Mat>.*),\s?Inc\s?:\s?(?P<Inc>.*)\s?,\s?rmean\s?:\s?(?P<rmean>[-+]?[0-9]*\.?[0-9]*)\s?,\s?phi\s?:\s?(?P<phi>[-+]?[0-9]*\.?[0-9]*)\s?,\s?poly\s?:\s?(?P<poly>.[-+]?[0-9]*\.?[0-9]*).*}|})'
    config_reader.add_pattern('layers', None, regex_layers, retrieve_layer)
    config_reader.add_pattern('halfspace_left', str, regex='default')
    config_reader.add_pattern('halfspace_right', str, regex='default')
    config_reader.add_pattern('plotphase', bool, regex='default')
    return config_reader.process_file()






def create_config_file(path_to_file, freq, angle, mat):
    """
    Function to write a simple config file template. For example:

        mat = ["eau"] + [("acier", 100)] + ["eau"]
        angle = [0]
        freq = np.linspace(0,1,101)

    """
    import os
    from datetime import date

    today = date.today().strftime("%d %B %Y")

    config_file = open(path_to_file,"w")

    # Small function to insert a EOL (end of line) character between each
    # line (each element of the "list of lines")
    ieol = lambda lst: [ lst[int(i/2)] if i%2==0 else "\n" for i in range(2*len(lst)) ]

    the_config = list()
    the_config.append("# Config file (" + today + ")")

    # definition of the layers
    for i, u in enumerate(mat[1:-1]):
        lay = "layer#{}=" . format(i+1) + "{"
        if u[0]=='meta':
            lay += "meta:{}" . format(u[1]) + ","
            lay += "Mat:{}" . format(u[2]) + ",Inc:{}" . format(u[3]) + ","
            lay += "rmean:{}" . format(u[4]/2) + ","
            lay += "phi:{}" . format(u[5]) + ",poly:{}" . format(u[6]) + "}"
        else:
            lay += u[0] + ":" + str(u[1]) + "}"
        the_config.append(lay)

    # definition of the left and right half spaces
    the_config.append("halfspace_left=" + mat[0])
    the_config.append("halfspace_right=" + mat[-1])

    # representation
    if len(freq)==1:
        the_config.append("theta_min={}" . format(angle[0]) )
        the_config.append("theta_max={}" . format(angle[-1]) )
        the_config.append("theta_num={}" . format(len(angle)) )
        the_config.append("f_fix  = {}" . format(freq[0]) )
        todo = 1
    elif len(angle)==1:
        the_config.append("f_min={}" . format(freq[0]) )
        the_config.append("f_max={}" . format(freq[-1]) )
        the_config.append("f_num={}" . format(len(freq)) )
        the_config.append("theta_fix={}" . format(angle[0]) )
        todo = 0
    else:
        pass

    the_config.append("todo=" + str(todo) )

    # write the config file with the list of line
    config_file.writelines( ieol(the_config) )
    config_file.close()

    return




if __name__ == '__main__':
    pass
