import {useSimulationStore} from '../state/simulationStore'
import {SimulationMode} from '../types'

const modes: { label: string; value: SimulationMode }[] = [
  { label: 'First Person', value: 'firstPerson' },
  { label: "Bird's-eye", value: 'birdseye' }
]

const ModeSwitcher = () => {
  const mode = useSimulationStore((s) => s.mode)
  const setMode = useSimulationStore((s) => s.setMode)
  return (
    <div className="button-row">
      {modes.map((m) => (
        <button
          key={m.value}
          className={`tab-button ${mode === m.value ? 'active' : ''}`}
          onClick={() => setMode(m.value)}
        >
          {m.label}
        </button>
      ))}
    </div>
  )
}

export default ModeSwitcher
