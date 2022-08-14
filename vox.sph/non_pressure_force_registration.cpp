//  Copyright (c) 2022 Feng Yang
//
//  I am making my contributions/submissions to this project solely in my
//  personal capacity and am not conveying any rights to any intellectual
//  property of any third parties.

#include "vox.sph/drag/drag_force_gissler2017.h"
#include "vox.sph/drag/drag_force_macklin2014.h"
#include "vox.sph/elasticity/elasticity_becker2009.h"
#include "vox.sph/elasticity/elasticity_kugelstadt2021.h"
#include "vox.sph/elasticity/elasticity_peer2018.h"
#include "vox.sph/simulation.h"
#include "vox.sph/surface_tension/surface_tension_akinci2013.h"
#include "vox.sph/surface_tension/surface_tension_becker2007.h"
#include "vox.sph/surface_tension/surface_tension_he2014.h"
#include "vox.sph/viscosity/viscosity_bender2017.h"
#include "vox.sph/viscosity/viscosity_peer2015.h"
#include "vox.sph/viscosity/viscosity_peer2016.h"
#include "vox.sph/viscosity/viscosity_standard.h"
#include "vox.sph/viscosity/viscosity_takahashi2015.h"
#include "vox.sph/viscosity/viscosity_weiler2018.h"
#include "vox.sph/viscosity/viscosity_xsph.h"
#include "vox.sph/vorticity/micropolar_model_bender2017.h"
#include "vox.sph/vorticity/vorticity_confinement.h"
#ifdef USE_THIRD_PARTY_METHODS
#include "SurfaceTension/SurfaceTension_ZorillaRitter2020.h"
#endif

using namespace vox;

void Simulation::registerNonpressureForces() {
    addDragMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
    addDragMethod("Macklin et al. 2014", DragForce_Macklin2014::creator);
    addDragMethod("Gissler et al. 2017", DragForce_Gissler2017::creator);

    addElasticityMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
    addElasticityMethod("Becker et al. 2009", Elasticity_Becker2009::creator);
    addElasticityMethod("Peer et al. 2018", Elasticity_Peer2018::creator);
    addElasticityMethod("Kugelstadt et al. 2021", Elasticity_Kugelstadt2021::creator);

    addSurfaceTensionMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
    addSurfaceTensionMethod("Becker & Teschner 2007", SurfaceTension_Becker2007::creator);
    addSurfaceTensionMethod("Akinci et al. 2013", SurfaceTension_Akinci2013::creator);
    addSurfaceTensionMethod("He et al. 2014", SurfaceTension_He2014::creator);
#ifdef USE_THIRD_PARTY_METHODS
    addSurfaceTensionMethod("Zorilla, Ritter, et al. 2020", SurfaceTension_ZorillaRitter2020::creator);
#endif

    addViscosityMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
    addViscosityMethod("Standard", Viscosity_Standard::creator);
    addViscosityMethod("XSPH", Viscosity_XSPH);
    addViscosityMethod("Bender and Koschier 2017", Viscosity_Bender2017::creator);
    addViscosityMethod("Peer et al. 2015", Viscosity_Peer2015::creator);
    addViscosityMethod("Peer et al. 2016", Viscosity_Peer2016::creator);
    addViscosityMethod("Takahashi et al. 2015 (improved)", Viscosity_Takahashi2015::creator);
    addViscosityMethod("Weiler et al. 2018", Viscosity_Weiler2018::creator);

    addVorticityMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
    addVorticityMethod("Micropolar model", MicropolarModel_Bender2017::creator);
    addVorticityMethod("Vorticity confinement", VorticityConfinement::creator);
}
