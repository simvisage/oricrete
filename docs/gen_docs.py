#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Dec 21, 2011 by: rch

from etsproxy.traits.api import \
    HasTraits, Str, Property, cached_property, \
    Enum

import os.path

from os.path import expanduser

HOME_DIR = expanduser("~")
# build directory
BUILD_DIR = os.path.join(HOME_DIR, '.oricrete', 'docs')
# output directory for the documentation
DOCS_DIR = os.path.join('..', 'docs',)
# output directory for the example documentation
EX_OUTPUT_DIR = os.path.join(DOCS_DIR, '_component_modules')

class GenExampleDoc(HasTraits):

    header = Str('''
%s
================================

    ''')

    component_module = None

    #===========================================================================
    # Derived traits
    #===========================================================================
    component_obj = Property(depends_on='component_module')
    @cached_property
    def _get_component_obj(self):
        return self.component_module.create_doc_object()

    name = Property(depends_on='component_module')
    @cached_property
    def _get_name(self):
        return self.component_obj.__class__

    output_dir = Property(depends_on='component_module')
    @cached_property
    def _get_output_dir(self):
        return os.path.join(EX_OUTPUT_DIR, self.name)

    rst_file_name = Property(depends_on='component_module')
    @cached_property
    def _get_rst_file_name(self):
        return os.path.join(self.output_dir, 'index.rst')

    def generate_html(self):

        print 'generating documentation for', self.name, '...'

        rst_text = '''
================================
Parametric study for %s
================================
        ''' % self.name

        dobj = self.component_obj

        if dobj.s.q.__doc__ != None:
            rst_text += dobj.s.q.__doc__

        rst_text += self.header

        for st in dobj.sampling_types:
            rst_text += '''

.. image:: %s_%s.png
    :width: 24%%

            ''' % (self.name, st)

        for st in dobj.sampling_types:
            rst_text += '''

.. image:: %s_sampling_%s.png
    :width: 24%%

            ''' % (self.name, st)

        rst_text += '\nFollowing oricrete configuration has been used to produce the sampling figures:\n\n'
        rst_text += '\n>>> print component_obj\n' + str(dobj.s) + '\n'

        rst_text += '''
Comparison of execution time for different sampling types
=========================================================
Execution time evaluated for an increasing number of sampling points n_sim:
'''
        for basename in dobj.fnames_sampling_efficiency:
            rst_text += '''

.. image:: %s
    :width: 100%%

            ''' % basename
            print 'written file %s', basename

        rst_text += '\n'

        rst_text += '''
Comparison of efficiency for different code types
=========================================================
Execution time evaluated for an numpy, weave and cython code:
'''
        for basename in dobj.fnames_language_efficiency:
            rst_text += '''

.. image:: %s
    :width: 100%%

            ''' % basename
            print 'written file %s', basename

        rst_text += '\n'

        rst_file = open(self.rst_file_name, 'w')

        rst_file.write(rst_text)

        rst_file.close()

class GenDoc(HasTraits):
    '''
    Configuration of the document generation using sphinx.
    '''
    component_modules = []

    build_mode = Enum('local', 'global')

    build_dir = Property(depends_on='build_mode')
    def _get_build_dir(self):
        build_dir = {'local' : '.',
                     'global' : BUILD_DIR }
        return build_dir[self.build_mode]

    html_server = 'root@mordred.imb.rwth-aachen.de:/var/www/docs/oricrete'

    method_dispatcher = {'all' : 'generate_examples' }

    def generate_html(self):
        for demo in self.component_modules:
            ged = GenExampleDoc(component_module=demo)
            ged.generate_html()

        os.chdir(DOCS_DIR)
        sphings_cmd = 'sphinx-build -b html -E . %s' % self.build_dir
        os.system(sphings_cmd)

    def push_html(self):
        '''
        Push the documentation to the server.
        '''
        rsync_cmd = 'rsync -av --delete %s/ %s' % (self.build_dir, self.html_server)
        os.system(rsync_cmd)

if __name__ == '__main__':

    gd = GenDoc(build_mode='global')
    #gd.generate_examples() # kind = 'sampling_efficiency')
    gd.generate_html()
    gd.push_html()
