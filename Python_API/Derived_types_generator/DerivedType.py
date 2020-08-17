class DerivedType():
    def __init__(self, python_type_name, fortran_type_name, fortran_module_name):
        self.python_type_name = python_type_name
        self.fortran_type_name = fortran_type_name
        self.fortran_module_name = fortran_module_name
        self.__params = []
    
    def addParam(self, param_python_name, param_fortran_name, type, description):
        # Types: string, object, enum, ndarray, class
        self.__params.append((param_python_name, param_fortran_name, type, description))
        
    def generate_template(self, file_name):
        self.__generate_python(file_name)
        self.__generate_init(file_name)
        self.__generate_fortran(file_name)
    
    def __generate_python(self, file_name):
        ret = ""
        for (param_python_name, param_fortran_name, type, description) in self.__params:
            if type == "bool":                
                ret +=\
    """
    def is_%s(self):""" % (param_python_name)
            else:
                ret +=\
    """
    @property
    def %s(self):"""  % (param_python_name)
    
            ret +=\
    """
        \"\"\"
        %s
        \"\"\"
        return crysfml_api.%s_get_%s(self.__address)["%s"]
    """ % (description, self.fortran_module_name, param_fortran_name, param_fortran_name)

        with(open(file_name, "a")) as file:
            file.write("========== Python part ==========")
            file.write(ret + "\n\n")
            
    def __generate_init(self, file_name):   
        ret = ""          
        for (_, param_fortran_name, type, _) in self.__params:
            ret += \
    """
    call method_table%%add_method("%s_get_%s", &                  ! method name
         "%s getter", &  !doc-string
         METH_VARARGS, &                  ! this method takes arguments but no keyword arguments
         c_funloc(%s_get_%s))  ! address of Fortran function to add
    """ % (self.fortran_module_name, param_fortran_name, param_fortran_name, self.fortran_module_name, param_fortran_name)
        
        with(open(file_name, "a")) as file:
            file.write("========== Fortran init part ==========")
            file.write(ret+"\n\n")
        
    def __generate_fortran(self, file_name):
        ret = ""
        ret += self.__generate_fortran_init(ret)
        ret += self.__generate_fortran_destructor(ret)
        ret += self.__generate_fortran_getters(ret)
        
        with(open(file_name, "a")) as file:
                    file.write("========== Fortran part ==========")
                    file.write(ret)
                    
    def __generate_fortran_init(self, ret):
        return ret
    
    def __generate_fortran_destructor(self, ret): 
        return ret
       
    def __generate_fortran_getters(self, ret):
        for (_, param_fortran_name, type, _) in self.__params:       
            ret += \
  """
  function %s_get_%s(self_ptr, args_ptr) result(r) bind(c)
        
    type(c_ptr), value :: self_ptr
    type(c_ptr), value :: args_ptr
    type(c_ptr) :: r
    type(tuple) :: args
    type(dict) :: retval
    integer :: num_args
    integer :: ierror
    type(%s_p) :: %s_pointer
  """ %(self.fortran_module_name, param_fortran_name, self.fortran_type_name, self.fortran_type_name.lower())
        
            if type == "ndarray":
                ret += \
    """
    type(ndarray) :: %s
    """ %(param_fortran_name)
        
            ret += \
    """
    r = C_NULL_PTR   ! in case of an exception return C_NULL_PTR
    ! use unsafe_cast_from_c_ptr to cast from c_ptr to tuple
    call unsafe_cast_from_c_ptr(args, args_ptr)
    ! Check if the arguments are OK
    ierror = args%%len(num_args)
    ! we should also check ierror, but this example does not do complete error checking for simplicity
    if (num_args /= 1) then
       call raise_exception(TypeError, "get_%s expects exactly 1 argument")
       call args%%destroy
       return
    endif
    
    ! Doing boring stuff
    call get_%s_from_arg(args, %s_pointer)""" % (param_fortran_name, self.fortran_type_name.lower(), self.fortran_type_name.lower())
            if type == "ndarray":
                ret += \
    """
    ierror = ndarray_create(%s, %s_pointer%%p%%%s)""" %  (param_fortran_name, self.fortran_type_name.lower(), param_fortran_name)
        
            ret += \
    """
    ierror = dict_create(retval)"""
            if type == "ndarray":
                ret += \
    """
    ierror = retval%%setitem("%s", %s)"""%  (param_fortran_name, param_fortran_name)
            elif type == "string":
                ret += \
    """
    ierror = retval%%setitem("%s", trim(%s_pointer%%p%%%s))
    """%  (param_fortran_name, self.fortran_type_name.lower(), param_fortran_name)        
            else:
                ret += \
    """
    ierror = retval%%setitem("%s", %s_pointer%%p%%%s)"""%  (param_fortran_name, self.fortran_type_name.lower(), param_fortran_name)
            
            ret += \
    """
    r = retval%%get_c_ptr()
    
  end function %s_get_%s
  
  """ %  (self.fortran_module_name, param_fortran_name)
  
        return ret
