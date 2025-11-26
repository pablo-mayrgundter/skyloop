import { useSimulationStore } from '../state/simulationStore'

const ControlsPanel = () => {
  const { isPlaying, simSpeed, play, pause, setSimSpeed, toggleOverlay } = useSimulationStore((s) => ({
    isPlaying: s.isPlaying,
    simSpeed: s.simSpeed,
    play: s.play,
    pause: s.pause,
    setSimSpeed: s.setSimSpeed,
    toggleOverlay: s.toggleOverlay
  }))

  const overlayButtons: { key: Parameters<typeof toggleOverlay>[0]; label: string }[] = [
    { key: 'showAtmosphere', label: 'Atmosphere' },
    { key: 'showWind', label: 'Wind' },
    { key: 'showTethers', label: 'Tethers' },
    { key: 'showSections', label: 'Sections' },
    { key: 'showGliders', label: 'Gliders' },
    { key: 'showDayNight', label: 'Day/Night' }
  ]

  return (
    <div className="panel">
      <h3>Controls</h3>
      <div className="button-row">
        <button onClick={isPlaying ? pause : play}>{isPlaying ? 'Pause' : 'Play'}</button>
        <button onClick={() => setSimSpeed(0.1)}>0.1x</button>
        <button onClick={() => setSimSpeed(1)}>1x</button>
        <button onClick={() => setSimSpeed(5)}>5x</button>
        <button onClick={() => setSimSpeed(10)}>10x</button>
      </div>
      <div style={{ marginTop: '8px' }}>Speed: {simSpeed.toFixed(2)}x</div>
      <h4>Overlays</h4>
      <div className="button-row">
        {overlayButtons.map((overlay) => (
          <button key={overlay.key} onClick={() => toggleOverlay(overlay.key)} className="tab-button">
            Toggle {overlay.label}
          </button>
        ))}
      </div>
    </div>
  )
}

export default ControlsPanel
