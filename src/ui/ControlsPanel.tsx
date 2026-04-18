import { useSimulationStore } from '../state/simulationStore'

const ControlsPanel = () => {
  const isPlaying = useSimulationStore((state) => state.isPlaying)
  const play = useSimulationStore((state) => state.play)
  const pause = useSimulationStore((state) => state.pause)
  const simSpeed = useSimulationStore((state) => state.simSpeed)
  const setSimSpeed = useSimulationStore((state) => state.setSimSpeed)
  const toggleOverlay = useSimulationStore((state) => state.toggleOverlay)

  return (
    <div className="controls-grid">
      <div className="button-row">
        <button className="tab-button" onClick={isPlaying ? pause : play}>
          {isPlaying ? 'Pause' : 'Play'}
        </button>
        <label style={{ display: 'grid', gap: '0.25rem' }}>
          Sim speed: {simSpeed.toFixed(1)}x
          <input
            type="range"
            min={0.1}
            max={5}
            step={0.1}
            value={simSpeed}
            onChange={(e) => setSimSpeed(Number(e.target.value))}
          />
        </label>
      </div>
      <div className="controls-grid">
        <div className="button-row">
          <button className="tab-button" onClick={() => toggleOverlay('showAtmosphere')}>
            Atmosphere
          </button>
          <button className="tab-button" onClick={() => toggleOverlay('showWind')}>
            Wind
          </button>
          <button className="tab-button" onClick={() => toggleOverlay('showDayNight')}>
            Day/Night
          </button>
        </div>
        <div className="button-row">
          <button className="tab-button" onClick={() => toggleOverlay('showSections')}>
            Sections
          </button>
          <button className="tab-button" onClick={() => toggleOverlay('showTethers')}>
            Tethers
          </button>
          <button className="tab-button" onClick={() => toggleOverlay('showGliders')}>
            Gliders
          </button>
        </div>
      </div>
    </div>
  )
}

export default ControlsPanel
