import { useSimulationStore } from '../state/simulationStore'

const SimulationHUD = () => {
  const { time, simSpeed, mode, isPlaying } = useSimulationStore((state) => ({
    time: state.time,
    simSpeed: state.simSpeed,
    mode: state.mode,
    isPlaying: state.isPlaying
  }))

  return (
    <div className="controls-grid">
      <div className="hud-row">
        <div className="metric">
          <div>Simulation Time</div>
          <strong>{time.toFixed(1)} s</strong>
        </div>
        <div className="metric">
          <div>Speed</div>
          <strong>{simSpeed.toFixed(1)}x</strong>
        </div>
        <div className="metric">
          <div>Mode</div>
          <strong>{mode === 'firstPerson' ? 'First Person' : "Bird's-eye"}</strong>
        </div>
        <div className="metric">
          <div>Status</div>
          <strong>{isPlaying ? 'Playing' : 'Paused'}</strong>
        </div>
      </div>
    </div>
  )
}

export default SimulationHUD
