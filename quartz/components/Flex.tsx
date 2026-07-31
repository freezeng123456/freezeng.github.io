import { concatenateResources } from "../util/resources"
import { classNames } from "../util/lang"
import { QuartzComponent, QuartzComponentConstructor, QuartzComponentProps } from "./types"
import LanguageSwitcher from "./LanguageSwitcher"

type FlexConfig = {
  components: {
    Component: QuartzComponent
    grow?: boolean
    shrink?: boolean
    basis?: string
    order?: number
    align?: "start" | "end" | "center" | "stretch"
    justify?: "start" | "end" | "center" | "between" | "around"
  }[]
  direction?: "row" | "row-reverse" | "column" | "column-reverse"
  wrap?: "nowrap" | "wrap" | "wrap-reverse"
  gap?: string
}

export default ((config: FlexConfig) => {
  const isToolbar =
    config.components.some((c) => String(c.Component.css).includes(".search")) &&
    config.components.some((c) => String(c.Component.css).includes(".darkmode"))
  const components = isToolbar
    ? [...config.components, { Component: LanguageSwitcher() }]
    : config.components

  const Flex: QuartzComponent = (props: QuartzComponentProps) => {
    const direction = config.direction ?? "row"
    const wrap = config.wrap ?? "nowrap"
    const gap = config.gap ?? "1rem"

    return (
      <div
        class={classNames(props.displayClass, "flex-component")}
        style={`flex-direction: ${direction}; flex-wrap: ${wrap}; gap: ${gap};`}
      >
        {components.map((c) => {
          const grow = c.grow ? 1 : 0
          const shrink = (c.shrink ?? true) ? 1 : 0
          const basis = c.basis ?? "auto"
          const order = c.order ?? 0
          const align = c.align ?? "center"
          const justify = c.justify ?? "center"

          return (
            <div
              style={`flex-grow: ${grow}; flex-shrink: ${shrink}; flex-basis: ${basis}; order: ${order}; align-self: ${align}; justify-self: ${justify};`}
            >
              <c.Component {...props} />
            </div>
          )
        })}
      </div>
    )
  }

  Flex.afterDOMLoaded = concatenateResources(...components.map((c) => c.Component.afterDOMLoaded))
  Flex.beforeDOMLoaded = concatenateResources(...components.map((c) => c.Component.beforeDOMLoaded))
  Flex.css = concatenateResources(...components.map((c) => c.Component.css))
  return Flex
}) satisfies QuartzComponentConstructor<FlexConfig>
