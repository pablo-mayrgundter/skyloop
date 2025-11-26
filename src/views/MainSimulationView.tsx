import MainCanvas from '../three/MainCanvas'
import ControlsPanel from '../ui/ControlsPanel'
import SimulationHUD from '../ui/SimulationHUD'
import ChartsPanel from '../ui/ChartsPanel'

const MainSimulationView = () => (
  <div className="main-content">
    <div className="panel">
      <div className="canvas-wrapper">
        <MainCanvas />
      </div>
    </div>
    <div className="panel" style={{ display: 'flex', flexDirection: 'column', gap: '12px' }}>
      <SimulationHUD />
      <ControlsPanel />
      <ChartsPanel />
    </div>
  </div>
)

export default MainSimulationView
