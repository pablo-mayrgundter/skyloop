import AtmosphereModule from '../modules/atmosphere/AtmosphereModule'
import BomModule from '../modules/bom/BomModule'
import DayNightModule from '../modules/dayNight/DayNightModule'
import EconomicsModule from '../modules/economics/EconomicsModule'
import GlidersModule from '../modules/gliders/GlidersModule'
import SectionsModule from '../modules/sections/SectionsModule'
import TethersModule from '../modules/tethers/TethersModule'
import TubeModule from '../modules/tube/TubeModule'
import WeatherModule from '../modules/weather/WeatherModule'
import {useSimulationStore} from '../state/simulationStore'
import ModuleSelector from '../ui/ModuleSelector'

const ModuleLabView = () => {
  const selected = useSimulationStore((s) => s.selectedModule)

  const moduleComponent = {
    sections: <SectionsModule />,
    tethers: <TethersModule />,
    tube: <TubeModule />,
    gliders: <GlidersModule />,
    atmosphere: <AtmosphereModule />,
    weather: <WeatherModule />,
    dayNight: <DayNightModule />,
    bom: <BomModule />,
    economics: <EconomicsModule />,
    main: <SectionsModule />
  }[selected]

  return (
    <div className="main-content">
      <div className="panel" style={{ width: '100%' }}>
        <div className="module-layout">
          <ModuleSelector />
          <div>{moduleComponent}</div>
        </div>
      </div>
    </div>
  )
}

export default ModuleLabView
