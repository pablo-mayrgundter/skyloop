import { useState } from 'react'
import ModuleLabView from './views/ModuleLabView'
import MainSimulationView from './views/MainSimulationView'
import ModeSwitcher from './ui/ModeSwitcher'

const App = () => {
  const [view, setView] = useState<'main' | 'lab'>('main')

  return (
    <div className="app-shell">
      <header className="header">
        <h1>Skyloop</h1>
        <button className={`tab-button ${view === 'main' ? 'active' : ''}`} onClick={() => setView('main')}>
          Main Simulation
        </button>
        <button className={`tab-button ${view === 'lab' ? 'active' : ''}`} onClick={() => setView('lab')}>
          Module Lab
        </button>
        <div style={{ marginLeft: 'auto' }}>
          <ModeSwitcher />
        </div>
      </header>
      {view === 'main' ? <MainSimulationView /> : <ModuleLabView />}
    </div>
  )
}

export default App
