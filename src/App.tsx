import { useState } from 'react'
import MainSimulationView from './views/MainSimulationView'
import ModuleLabView from './views/ModuleLabView'
import ModeSwitcher from './ui/ModeSwitcher'
import { SimulationMode } from './types/simulation'
import { useSimulationStore } from './state/simulationStore'

const App = () => {
  const [activeView, setActiveView] = useState<'main' | 'lab'>('main')
  const mode = useSimulationStore((state) => state.mode)
  const setMode = useSimulationStore((state) => state.setMode)

  const renderView = () => {
    if (activeView === 'lab') {
      return <ModuleLabView />
    }
    return <MainSimulationView />
  }

  return (
    <div className="app-shell">
      <header className="header">
        <div className="brand">
          <span>Skyloop</span>
          <span style={{ opacity: 0.7, fontSize: '0.9rem' }}>Earth-to-Earth Transit</span>
        </div>
        <div className="button-row">
          <button
            className={`tab-button ${activeView === 'main' ? 'active' : ''}`}
            onClick={() => setActiveView('main')}
          >
            Main Simulation
          </button>
          <button
            className={`tab-button ${activeView === 'lab' ? 'active' : ''}`}
            onClick={() => setActiveView('lab')}
          >
            Module Lab
          </button>
          <ModeSwitcher
            value={mode}
            onChange={(newMode: SimulationMode) => setMode(newMode)}
          />
        </div>
      </header>
      {renderView()}
    </div>
  )
}

export default App
