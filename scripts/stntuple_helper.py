#!/usr/bin/env python
#
#
import os, re, string, subprocess
#------------------------------------------------------------------------------
class stntuple_helper:
    """mu2e_helper: class to produce library names"""
#   This appears to behave like c++ static member and is initialized at class defintion time.
    sourceroot =  os.path.abspath('.')

    def __init__(self, env, debug=False):
        self._env    = env;
        self._list_of_object_files = [];

        self.dd      = re.search('[^/]*/[^/]*$',self._env.Dir('.').abspath).group(0)
        self.dirname = os.path.dirname(self.dd);   # THIS
        self.subdir  = os.path.basename(self.dd);
        self.libname = self.dirname+'_'+self.subdir
        self.d1      = self.libname+'-shared';
        self.suffix  = ".hh" ;
        self.tmp_dir = "tmp/src/"+self.d1;
        self._debug  = debug
        if (debug) : print "-------------- building directory: "+self.dirname+'/'+self.subdir
#
#   Accesor
#
    def base(self):
        return self.sourceroot

    def tmp_dira(self):
        return self.tmp_dir

    def debug(self):
        return self._debug

    def do_codegen(self,cc,script):
        self._env.StntupleCodegen(cc,script);
        obj = string.replace(cc,".cc",".o");
        self._env.SharedObject(obj,cc)
        self._list_of_object_files.append(obj);

    def handle_dictionaries(self,suffix = ".hh",skip_list = []):
        self.suffix = suffix ;
#------------------------------------------------------------------------------
# generate dictionaries
#------------------------------------------------------------------------------
        list_of_linkdef_files = self._env.Glob(self.subdir+'/dict/*_linkdef.h', strings=True)
        if (self._debug):
            print ("["+self.dirname+"/"+self.subdir+"] handle_dictionaries: list_of_linkdef_files = ",list_of_linkdef_files)
            
        list_of_dict_files    = []

        for f in list_of_linkdef_files:
            linkdef       = string.split(str(f),'/');
            linkdef_fn    = linkdef[len(linkdef)-1];

            if (self._debug): 
                print ("linkdef_fn = ",linkdef_fn)

            if (not linkdef_fn in skip_list):
                clname        = string.replace(linkdef_fn,"_linkdef.h","");
                include       = self.subdir+'/'+clname+self.suffix;
            
                dict          = '#/tmp/src/'+self.d1+'/'+clname+'_dict.cxx';
                list_of_dict_files.append(dict);

                if (self._debug):
                    print ("linkdef = ",linkdef)
                    print ("include = ",include)

                self._env.StntupleRootCint(dict,[f,include])
#------------------------------------------------------------------------------
# compile dictionaries
#------------------------------------------------------------------------------
        list   = [];

        for dict in list_of_dict_files:
            obj_cxx_file = string.replace(dict,".cxx",".o");
            list.append(obj_cxx_file);
            self._list_of_object_files.append(obj_cxx_file);
        
            include  = string.replace(dict,".cxx",".h");
            self._env.SharedObject(obj_cxx_file,dict)


    def compile_fortran(self,list_of_f_files, skip_list = []):
        if (self._debug):
            print  (self._dirname+"[build_libs]: list_of_f_files:"+self._subdir,list_of_f_files);

        for f in list_of_f_files:
            if (not f in skip_list):
                o = '#tmp/src/'+self.d1+'/'+string.split(f,'.')[0]+'.o'
                self._list_of_object_files.append(o);
                self._env.SharedObject(o,f)


    def build_libs(self,list_of_cc_files, skip_list = [],libs = []):
        if (self._debug):
            print ("["+self.dirname+".build_libs]: list_of_cc_files:"+self.subdir,list_of_cc_files)

        for cc in list_of_cc_files:
            if (not cc in skip_list):
                if (self._debug):
                    print (".cc file: "+cc)
                    
                o = '#/tmp/src/'+self.d1+'/'+string.split(cc,'.')[0]+'.o'
                self._list_of_object_files.append(o);
                self._env.SharedObject(o,cc)

        dir = os.environ.get('MU2E_SATELLITE_RELEASE');
        if (dir == None) : dir = os.environ.get('MU2E_BASE_RELEASE');
        lib_name = dir+'/lib/'+self.libname+'.so';

        if (self._debug):
            print ("list_of_object_files:",self._list_of_object_files)

        self._env.SharedLibrary(lib_name,self._list_of_object_files,LIBS = [libs])


    def build_modules(self,list_of_module_files, skip_list, libs = []):
        if (self._debug):
            print ("[Stntuple.build_libs]: list_of_module_files in Stntuple/"+self.subdir," : ",list_of_module_files)

        for module in list_of_module_files:
            if (not module in skip_list):
                if (self._debug):
                    print ("module file: "+module)
                    
                o = '#/tmp/src/'+self.d1+'/'+string.split(module,'.')[0]+'.o'
                self._env.SharedObject(o,module)

                mname = string.split(os.path.basename(module),'.')[0];
                lib   = '#/lib/libmu2e_'+self.dirname+'_'+mname+'.so';
                if (self._debug):
                    print ("o: "+o, "lib:"+lib)
                    
                self._env.SharedLibrary(lib,o,LIBS =  ['libStntuple_mod.so',libs]);
#------------------------------------------------------------------------------ END
