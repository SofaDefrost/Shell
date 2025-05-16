#ifndef SOFA_GUI_ADAPTIVECUTTING_H
#define SOFA_GUI_ADAPTIVECUTTING_H

#include <sofa/component/setting/MouseButtonSetting.h>
#include <sofa/gui/component/performer/MouseInteractor.h>
#include <sofa/gui/component/performer/InteractionPerformer.h>
#include <sofa/gui/common/MouseOperations.h>
#include <sofa/gui/common/OperationFactory.h>

#include <Shell/Adaptivity/controller/AdaptiveCuttingController.h>

namespace sofa
{

namespace component
{

// ----------  Settings  -----------------------------------------------------

namespace configurationsetting
{

class SHELL_ADAPTIVITY_API AdaptiveCuttingSetting: public sofa::component::setting::MouseButtonSetting
{
public:
    SOFA_CLASS(AdaptiveCuttingSetting, sofa::component::setting::MouseButtonSetting);

protected:
    AdaptiveCuttingSetting() { }

public:
    std::string getOperationType() override { return "AdaptiveCutting"; }
};

} // namespace configurationsetting

// ----------  Performer  ----------------------------------------------------

namespace collision
{

template <class DataTypes>
class AdaptiveCuttingPerformer: public sofa::gui::component::performer::TInteractionPerformer<DataTypes>
{
public:
    //typedef sofa::component::collision::BaseContactMapper< DataTypes >        MouseContactMapper;
    //typedef sofa::core::behavior::MechanicalState< DataTypes >         MouseContainer;
    //typedef sofa::core::behavior::BaseForceField              MouseForceField;

    AdaptiveCuttingPerformer(gui::component::performer::BaseMouseInteractor *i) : sofa::gui::component::performer::TInteractionPerformer<DataTypes>(i)
    {}

    void start() {
        //std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl;
        /*TODO*/
    }

    void execute() {
        //std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl;
        //updatePosition();
        /* TODO */
    }

    //void updatePosition() {
    //    if (!this->interactor) return;
    //    BodyPicked picked = this->interactor->getBodyPicked();

    //    if (!picked.body && !picked.mstate) return;
    //    core::behavior::MechanicalState<DataTypes>* mstate=NULL;

    //    if (picked.body) {
    //        mstate = dynamic_cast<core::behavior::MechanicalState<DataTypes>*>(
    //           picked.body->getContext()->getMechanicalState());
    //    } else {
    //        mstate = dynamic_cast<core::behavior::MechanicalState<DataTypes>*>(
    //            picked.mstate);
    //    }

    //    if (!mstate)
    //    {
    //        this->interactor->serr
    //            << "Incompatible MState during Mouse Interaction"
    //            << this->interactor->sendl;
    //        return;
    //    }

    //    sofa::component::controller::Test2DAdapter<DataTypes>* adapter;
    //    mstate->getContext()->get(adapter);

    //    if (!adapter) return;
    //    adapter->setTrackedPoint(picked);
    //}


    void configure(sofa::component::setting::MouseButtonSetting* setting)
    {
        configurationsetting::AdaptiveCuttingSetting* s = dynamic_cast<configurationsetting::AdaptiveCuttingSetting*>(setting);
        if (s)
        {
             //std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl; /*TODO*/
        }
    }
};

} // namespace collision

} // namespace component


// ----------  Operation  ----------------------------------------------------

namespace gui
{

class SHELL_ADAPTIVITY_API AdaptiveCuttingOperation : public sofa::gui::common::Operation
{
public:

    typedef ::sofa::component::controller::CuttingAdapter CuttingAdapter;

    AdaptiveCuttingOperation(sofa::component::configurationsetting::AdaptiveCuttingSetting::SPtr s = sofa::core::objectmodel::New<sofa::component::configurationsetting::AdaptiveCuttingSetting>()) : Operation(s)
    {}
    //virtual ~AdaptiveCuttingOperation() {}

    /// This function is called each time the mouse is clicked.
    void start();
    //void execution() { std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl; /* TODO */ Operation::execution();}
    /// This function is called after each mouse click.
    //void end(); { std::cout << __PRETTY_FUNCTION__ << " not implemented" << std::endl; /* TODO */ Operation::end();}
    /// This function is called when shift key is released.
    void endOperation();
    /// Called on mouse move
    void wait();

    static std::string getDescription() {return "Cutting with adaptive mesh";}

    std::string defaultPerformerType() { return "AdaptiveCutting"; }

private:

    /// Get the adaptive component associated with selected object.
    CuttingAdapter* getAdapter();


};

} // namespace gui

} // namespace sofa

#endif // #ifndef SOFA_GUI_ADAPTIVECUTTING_H
