import { Canvas } from '@react-three/fiber'
import { OrbitControls } from '@react-three/drei'
import LoopTethers from '../../three/loop/LoopTethers'

const TethersModule = () => {
  return (
    <div className="canvas-wrapper" style={{ height: 320 }}>
      <Canvas camera={{ position: [0, 0.5, 1] }}>
        <ambientLight intensity={0.8} />
        <LoopTethers />
        <OrbitControls />
      </Canvas>
      <p>Analyzing tether tension, damping, and ground coupling strategies.</p>
    </div>
  )
}

export default TethersModule
