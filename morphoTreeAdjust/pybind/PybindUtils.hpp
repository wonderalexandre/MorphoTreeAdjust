
#include "../include/Common.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#ifndef PYBIND_UTILS_H
#define PYBIND_UTILS_H

class PybindUtils{
    public:
        
        template <typename PixelType>
        static py::array_t<PixelType> toNumpy(ImagePtr<PixelType> image) {
            
            std::shared_ptr<PixelType[]> buffer = image->rawDataPtr();
            int n = image->getSize();
            std::shared_ptr<PixelType[]> bufferCopy = buffer;

            py::capsule free_when_done(new std::shared_ptr<PixelType[]>(bufferCopy), [](void* ptr) {
                // Converte de volta e destrói corretamente
                delete reinterpret_cast<std::shared_ptr<PixelType[]>*>(ptr);
            });
            
            py::array_t<PixelType> numpy = py::array(py::buffer_info(
                buffer.get(),
                sizeof(PixelType),
                py::format_descriptor<PixelType>::value,
                1,
                { n },
                { sizeof(PixelType) }
            ), free_when_done);
            
            return numpy;
        }
    


};

#endif
