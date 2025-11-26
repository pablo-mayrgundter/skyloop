import { useSimulationStore } from '../state/simulationStore'
import { SkyloopModuleId } from '../types/simulation'

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
  const selectedModule = useSimulationStore((state) => state.selectedModule)
  const setSelectedModule = useSimulationStore((state) => state.setSelectedModule)

  return (
    <div className="module-list">
      {modules.map((module) => (
        <button
          key={module.id}
          className={`tab-button module-item ${selectedModule === module.id ? 'active' : ''}`}
          onClick={() => setSelectedModule(module.id)}
        >
          {module.label}
        </button>
      ))}
    </div>
  )
}

export default ModuleSelector
