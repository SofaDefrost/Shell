#include <Shell/Adaptivity/init.h>
#include <sofa/core/ObjectFactory.h>

namespace shelladaptivity
{

void initializePlugin() 
{
    static bool first = true;
    if (first) {
        first = false;
        // Register components here
    }
}

}

extern "C" 
{
    SHELL_ADAPTIVITY_API void initExternalModule() 
    {
        shelladaptivity::initializePlugin();
    }

    SHELL_ADAPTIVITY_API const char* getModuleName() 
    {
        return shelladaptivity::MODULE_NAME;
    }

    SHELL_ADAPTIVITY_API const char* getModuleVersion() 
    {
        return shelladaptivity::MODULE_VERSION;
    }

    SHELL_ADAPTIVITY_API const char* getModuleLicense() 
    {
        return "";
    }

    SHELL_ADAPTIVITY_API const char* getModuleDescription() 
    {
        return "SOFA plugin for Shell.Adaptivity";
    }
}
