import MainCanvas from '../three/MainCanvas'
import SimulationHUD from '../ui/SimulationHUD'
import ControlsPanel from '../ui/ControlsPanel'
import ChartsPanel from '../ui/ChartsPanel'

const MainSimulationView = () => {
  return (
    <div className="main-grid">
      <div className="canvas-wrapper">
        <MainCanvas />
      </div>
      <div className="panel">
        <SimulationHUD />
        <ControlsPanel />
        <ChartsPanel />
      </div>
    </div>
  )
}

export default MainSimulationView
