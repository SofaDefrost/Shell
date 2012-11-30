#ifndef SOFA_GUI_ADAPTIVECUTTING_H
#define SOFA_GUI_ADAPTIVECUTTING_H

//#include <sofa/core/objectmodel/ConfigurationSetting.h>
#include <sofa/component/component.h>
#include <sofa/component/configurationsetting/MouseButtonSetting.h>
#include <sofa/component/collision/MouseInteractor.h>
//#include <sofa/component/collision/BaseContactMapper.h>
#include <sofa/component/collision/InteractionPerformer.h>
#include <sofa/gui/MouseOperations.h>

namespace sofa
{

namespace component
{

// ----------  Settings  -----------------------------------------------------

namespace configurationsetting
{

class SOFA_SHELLS_API AdaptiveCuttingSetting: public MouseButtonSetting
{
public:
    SOFA_CLASS(AdaptiveCuttingSetting, MouseButtonSetting);

protected:
    AdaptiveCuttingSetting() { }

public:
    std::string getOperationType() { return "AdaptiveCutting"; }
};

} // namespace configurationsetting

// ----------  Performer  ----------------------------------------------------

namespace collision
{

template <class DataTypes>
class AdaptiveCuttingPerformer: public TInteractionPerformer<DataTypes>
{
public:
    //typedef sofa::component::collision::BaseContactMapper< DataTypes >        MouseContactMapper;
    //typedef sofa::core::behavior::MechanicalState< DataTypes >         MouseContainer;
    //typedef sofa::core::behavior::BaseForceField              MouseForceField;

    AdaptiveCuttingPerformer(BaseMouseInteractor *i) : TInteractionPerformer<DataTypes>(i)
    {}

    void start() { std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl; /*TODO*/}
    void execute() { std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl; /*TODO*/}


    void configure(configurationsetting::MouseButtonSetting* setting)
    {
        configurationsetting::AdaptiveCuttingSetting* s = dynamic_cast<configurationsetting::AdaptiveCuttingSetting*>(setting);
        if (s)
        {
             std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl; /*TODO*/
        }
    }
};

} // namespace collision

} // namespace component


// ----------  Operation  ----------------------------------------------------

namespace gui
{

class SOFA_SHELLS_API AdaptiveCuttingOperation : public Operation
{
public:

    AdaptiveCuttingOperation(sofa::component::configurationsetting::AdaptiveCuttingSetting::SPtr s = sofa::core::objectmodel::New<sofa::component::configurationsetting::AdaptiveCuttingSetting>()) : Operation(s)
    {}
    //virtual ~AdaptiveCuttingOperation() {}

    /// This function is called each time the mouse is clicked.
    void start() { std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl; /* TODO */ Operation::start(); }
    void execution() { std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl; /* TODO */ Operation::execution();}
    /// This function is called after each mouse click.
    void end() { std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl; /* TODO */ Operation::end();}
    /// This function is called when shift key is released.
    void endOperation() {  { std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl; /* TODO */}; this->end(); }
    void wait() { { std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl; /* TODO */}}

    static std::string getDescription() {return "Cutting with adaptive mesh";}

    std::string defaultPerformerType() { return "AdaptiveCutting"; }

private:

};

} // namespace gui

} // namespace sofa

#endif // #ifndef SOFA_GUI_ADAPTIVECUTTING_H
