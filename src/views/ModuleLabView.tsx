import ModuleSelector from '../ui/ModuleSelector'
import SectionsModule from '../modules/sections/SectionsModule'
import TethersModule from '../modules/tethers/TethersModule'
import TubeModule from '../modules/tube/TubeModule'
import GlidersModule from '../modules/gliders/GlidersModule'
import AtmosphereModule from '../modules/atmosphere/AtmosphereModule'
import WeatherModule from '../modules/weather/WeatherModule'
import DayNightModule from '../modules/dayNight/DayNightModule'
import BomModule from '../modules/bom/BomModule'
import EconomicsModule from '../modules/economics/EconomicsModule'
import { useSimulationStore } from '../state/simulationStore'

const ModuleLabView = () => {
  const selectedModule = useSimulationStore((state) => state.selectedModule)

  const renderModule = () => {
    switch (selectedModule) {
      case 'sections':
        return <SectionsModule />
      case 'tethers':
        return <TethersModule />
      case 'tube':
        return <TubeModule />
      case 'gliders':
        return <GlidersModule />
      case 'atmosphere':
        return <AtmosphereModule />
      case 'weather':
        return <WeatherModule />
      case 'dayNight':
        return <DayNightModule />
      case 'bom':
        return <BomModule />
      case 'economics':
        return <EconomicsModule />
      default:
        return <SectionsModule />
    }
  }

  return (
    <div className="main-grid">
      <div className="panel">
        <h3>Modules</h3>
        <ModuleSelector />
      </div>
      <div className="panel">{renderModule()}</div>
    </div>
  )
}

export default ModuleLabView
