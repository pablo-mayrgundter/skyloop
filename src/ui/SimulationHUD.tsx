import { useMemo } from 'react'
import { SimulationMetrics } from '../types'
import { useSimulationStore } from '../state/simulationStore'

const formatTime = (seconds: number) => seconds.toFixed(1)

const SimulationHUD = () => {
  const { mode, simSpeed, time } = useSimulationStore((s) => ({
    mode: s.mode,
    simSpeed: s.simSpeed,
    time: s.time
  }))

  const metrics: SimulationMetrics = useMemo(() => ({
    velocityMps: 300,
    altitudeKm: 30,
    tetherTensionKn: 1200,
    powerMw: 45
  }), [])

  return (
    <div className="panel">
      <h3>Simulation HUD</h3>
      <div className="hud-row">
        <div className="metric-card">
          <strong>Mode</strong>
          <div>{mode}</div>
        </div>
        <div className="metric-card">
          <strong>Sim Speed</strong>
          <div>{simSpeed.toFixed(2)}x</div>
        </div>
        <div className="metric-card">
          <strong>Time</strong>
          <div>{formatTime(time)} s</div>
        </div>
      </div>
      <div className="hud-row" style={{ marginTop: '8px' }}>
        <div className="metric-card">
          <strong>Velocity</strong>
          <div>{metrics.velocityMps} m/s</div>
        </div>
        <div className="metric-card">
          <strong>Altitude</strong>
          <div>{metrics.altitudeKm} km</div>
        </div>
        <div className="metric-card">
          <strong>Tether Tension</strong>
          <div>{metrics.tetherTensionKn} kN</div>
        </div>
        <div className="metric-card">
          <strong>Power</strong>
          <div>{metrics.powerMw} MW</div>
        </div>
      </div>
    </div>
  )
}

export default SimulationHUD
