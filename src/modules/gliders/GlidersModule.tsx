import { Canvas } from '@react-three/fiber'
import { OrbitControls } from '@react-three/drei'
import LoopGliders from '../../three/loop/LoopGliders'
import { useSimulationStore } from '../../state/simulationStore'

const GlidersModule = () => {
  const play = useSimulationStore((state) => state.play)
  const pause = useSimulationStore((state) => state.pause)
  const isPlaying = useSimulationStore((state) => state.isPlaying)

  return (
    <div className="canvas-wrapper" style={{ height: 320 }}>
      <Canvas camera={{ position: [0, 0.25, 0.7] }}>
        <ambientLight intensity={0.8} />
        <LoopGliders />
        <OrbitControls />
      </Canvas>
      <div className="button-row">
        <button className="tab-button" onClick={isPlaying ? pause : play}>
          {isPlaying ? 'Pause glider' : 'Play glider'}
        </button>
        <p>Stubbed glider physics will evolve into full rapier constraints.</p>
      </div>
    </div>
  )
}

export default GlidersModule
