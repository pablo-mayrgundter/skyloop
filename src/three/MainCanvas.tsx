import {Physics} from '@react-three/rapier'
import {Canvas, useFrame} from '@react-three/fiber'
import {Suspense} from 'react'
import {Color} from 'three'
import {useSimulationStore} from '../state/simulationStore'
import SceneRoot from './SceneRoot'

const SimulationTicker = () => {
  const {isPlaying, setTime, simSpeed, time} = useSimulationStore((s) => ({
    isPlaying: s.isPlaying,
    setTime: s.setTime,
    simSpeed: s.simSpeed,
    time: s.time
  }))
  useFrame((_, delta) => {
    if (isPlaying) {
      setTime(time + delta * simSpeed)
    }
  })
  return null
}

const MainCanvas = () => (
  <Canvas shadows gl={{ antialias: true }} onCreated={({ scene }) => scene.background = new Color('#030711')}>
    <color attach="background" args={[0x030711]} />
    <ambientLight intensity={0.3} />
    <directionalLight position={[1, 1, 0.5]} intensity={1} castShadow />
    <Suspense fallback={null}>
      <Physics>
        <SimulationTicker />
        <SceneRoot />
      </Physics>
    </Suspense>
  </Canvas>
)

export default MainCanvas
