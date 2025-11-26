import { useSimulationStore } from '../state/simulationStore'
import { SkyloopModuleId } from '../types'

const modules: { id: SkyloopModuleId; label: string }[] = [
  { id: 'sections', label: 'Sections' },
  { id: 'tethers', label: 'Tethers' },
  { id: 'tube', label: 'Tube' },
  { id: 'gliders', label: 'Gliders' },
  { id: 'atmosphere', label: 'Atmosphere' },
  { id: 'weather', label: 'Weather' },
  { id: 'dayNight', label: 'Day/Night' },
  { id: 'bom', label: 'Bill of Materials' },
  { id: 'economics', label: 'Economics' }
]

const ModuleSelector = () => {
  const selected = useSimulationStore((s) => s.selectedModule)
  const setSelectedModule = useSimulationStore((s) => s.setSelectedModule)
  return (
    <ul className="module-list">
      {modules.map((module) => (
        <li key={module.id}>
          <button
            className={`tab-button ${selected === module.id ? 'active' : ''}`}
            onClick={() => setSelectedModule(module.id)}
          >
            {module.label}
          </button>
        </li>
      ))}
    </ul>
  )
}

export default ModuleSelector
