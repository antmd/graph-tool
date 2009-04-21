bool __exception_thrown = false;
string __exception_error;

//variable definitions
${var_defs}

try
{
    //variable extraction
    ${var_extract}

    // where the actual code is included
    ${code}
}
catch (const GraphException& e)
{
    __exception_error = e.what();
    __exception_thrown = true;
}
catch (const bad_any_cast& e)
{
    __exception_error = e.what();
    __exception_error += " (wrong property map type?)";
    __exception_thrown = true;
}
catch (const std::exception& e)
{
    __exception_error = "unknown exception thrown: ";
    __exception_error += e.what();
    __exception_thrown = true;
}

python::dict return_vals;
return_vals["__exception_error"] = __exception_error;
return_vals["__exception_thrown"] = __exception_thrown;

// updated values will be inserted in return_vals below
${return_vals}

return_val = py::object(return_vals.ptr());
