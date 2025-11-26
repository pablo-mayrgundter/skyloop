import { SimulationMode } from '../types/simulation'

const ModeSwitcher = ({ value, onChange }: { value: SimulationMode; onChange: (mode: SimulationMode) => void }) => {
  return (
    <div className="button-row">
      <button
        className={`tab-button ${value === 'firstPerson' ? 'active' : ''}`}
        onClick={() => onChange('firstPerson')}
      >
        First Person
      </button>
      <button
        className={`tab-button ${value === 'birdseye' ? 'active' : ''}`}
        onClick={() => onChange('birdseye')}
      >
        Bird's-eye
      </button>
    </div>
  )
}

export default ModeSwitcher
