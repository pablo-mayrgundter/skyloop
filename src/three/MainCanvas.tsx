import { Canvas } from '@react-three/fiber'
import { OrbitControls, StatsGl } from '@react-three/drei'
import SceneRoot from './SceneRoot'
import { PhysicsProvider } from '../physics'
import { useSimulationStore } from '../state/simulationStore'

const MainCanvas = () => {
  const mode = useSimulationStore((state) => state.mode)

  return (
    <Canvas shadows camera={{ position: [0, 1.2, 3], fov: 55 }}>
      <color attach="background" args={[0.02, 0.04, 0.09]} />
      <ambientLight intensity={0.2} />
      <directionalLight position={[3, 2, 1]} intensity={1.2} />
      <PhysicsProvider mode={mode}>
        <SceneRoot />
      </PhysicsProvider>
      <OrbitControls enablePan={false} />
      <StatsGl />
    </Canvas>
  )
}

export default MainCanvas
