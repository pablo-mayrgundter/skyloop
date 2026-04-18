import { Canvas } from '@react-three/fiber'
import LoopSections from '../../three/loop/LoopSections'
import { OrbitControls } from '@react-three/drei'

const SectionsModule = () => {
  return (
    <div className="canvas-wrapper" style={{ height: 320 }}>
      <Canvas camera={{ position: [0, 0.3, 0.6] }}>
        <ambientLight intensity={0.5} />
        <LoopSections />
        <OrbitControls />
      </Canvas>
      <p>Sections module focuses on buoyant vacuum segments along the great-circle arc.</p>
    </div>
  )
}

export default SectionsModule
