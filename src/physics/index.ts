import { ReactNode, useEffect } from 'react'
import { RigidBody, World, init, RigidBodyDesc } from '@react-three/rapier'
import { useFrame } from '@react-three/fiber'
import { useSimulationStore } from '../state/simulationStore'

let rapierWorld: World | null = null

const ensureWorld = async () => {
  if (rapierWorld) return rapierWorld
  await init()
  rapierWorld = new World({ x: 0, y: -9.81, z: 0 })
  return rapierWorld
}

export const usePhysicsBody = (descriptor: RigidBodyDesc) => {
  useEffect(() => {
    ensureWorld()
    return () => {}
  }, [descriptor])
  return { ref: undefined as unknown as React.RefObject<RigidBody> }
}

export const PhysicsProvider = ({ children, mode }: { children: ReactNode; mode: string }) => {
  const { isPlaying, simSpeed, time, setTime } = useSimulationStore((state) => ({
    isPlaying: state.isPlaying,
    simSpeed: state.simSpeed,
    time: state.time,
    setTime: state.setTime
  }))

  useFrame((_, delta) => {
    if (!rapierWorld) return
    if (isPlaying) {
      rapierWorld.timestep = delta * simSpeed
      rapierWorld.step()
      setTime(time + delta * simSpeed)
    }
  })

  return <>{children}</>
}

export const createGliderBody = async () => {
  const world = await ensureWorld()
  const bodyDesc = RigidBodyDesc.dynamic().setTranslation(0, 0, 0)
  const body = world.createRigidBody(bodyDesc)
  // TODO: attach collider and aerodynamic model
  return body
}

export const createSectionBody = async (position: { x: number; y: number; z: number }) => {
  const world = await ensureWorld()
  const bodyDesc = RigidBodyDesc.fixed().setTranslation(position.x, position.y, position.z)
  return world.createRigidBody(bodyDesc)
}

export const stepPhysics = async (deltaTime: number) => {
  const world = await ensureWorld()
  world.timestep = deltaTime
  world.step()
}
