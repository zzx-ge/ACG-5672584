#include <OpenImageDenoise/oidn.hpp>
#include "Imaging.h"
#include <iostream>

void denoiseImage(Colour* hdrpixels, int width, int height) {
    oidn::DeviceRef device = oidn::newDevice();
    device.commit();

    // Create the OIDN filter
    oidn::FilterRef filter = device.newFilter("RT"); // "RT" mode for ray tracing
    filter.setImage("color", hdrpixels, oidn::Format::Float3, width, height);
    filter.setImage("output", hdrpixels, oidn::Format::Float3, width, height);
    filter.set("hdr", true); // HDR images require this
    filter.commit();

    // Run the denoiser
    filter.execute();

    // Error handling
    const char* errorMessage;
    if (device.getError(errorMessage) != oidn::Error::None)
        std::cerr << "OIDN Error: " << errorMessage << std::endl;
}