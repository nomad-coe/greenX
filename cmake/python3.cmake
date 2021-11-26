
find_package(Python3 3.7 COMPONENTS Interpreter Development)
if(Python3_FOUND)
    message("-- Python 3 interpreter version: " ${Python3_VERSION})
else()
    message("-- Python 3 interpreter not found")
endif()
